"""BeamDeflection - Streamlit dashboard.

Run:  streamlit run app.py
Single-span beams (4 support types) and continuous beams (any number of spans, overhangs, clamped ends).
"""
import io
import math

import numpy as np
import pandas as pd
import streamlit as st

import beam_core as bc
import beam_plots as bp

st.set_page_config(page_title="BeamDeflection", layout="wide")
st.title("BeamDeflection dashboard")
st.caption("Linear-elastic Euler-Bernoulli beam, constant section. Not a substitute for a code-compliant design.")

# ----------------------------------------------------------------- sidebar: model
sb = st.sidebar
mode = sb.radio("Beam type", ["Single span", "Continuous beam"])

sb.subheader("Material")
mat_names = [m.name for m in bc.MATERIALS] + ["Custom"]
mat_name = sb.selectbox("Material", mat_names, index=1)
if mat_name == "Custom":
    e_gpa = sb.number_input("E (GPa)", 1.0, 500.0, 210.0)
    density = sb.number_input("Density (kg/m3)", 100.0, 20000.0, 7850.0)
    strength = sb.number_input("Strength (MPa)", 0.0, 2000.0, 275.0)
    shear_strength = 0.6 * strength if strength > 0 else None
    strength = strength if strength > 0 else None
else:
    mat = next(m for m in bc.MATERIALS if m.name == mat_name)
    e_gpa, density, strength, shear_strength = mat.e_gpa, mat.density, mat.strength_mpa, mat.shear_strength_mpa
    sb.write(f"E = {e_gpa} GPa, strength = {strength if strength else 'n/a'} MPa")

sb.subheader("Section")
sec_kind = sb.selectbox("Section type", ["IPE catalogue", "Rectangle", "Round bar", "Round tube", "Custom I and c"])
mass_per_m = None
if sec_kind == "IPE catalogue":
    sec_name = sb.selectbox("IPE size", [r[0] for r in bc.IPE], index=5)
    sec = bc.ipe(sec_name)
    mass_per_m = sec.mass_kg_per_m
elif sec_kind == "Rectangle":
    sec = bc.rect(sb.number_input("Width b (mm)", 5.0, 2000.0, 100.0), sb.number_input("Height h (mm)", 5.0, 3000.0, 200.0))
elif sec_kind == "Round bar":
    sec = bc.round_bar(sb.number_input("Diameter d (mm)", 5.0, 1000.0, 60.0))
elif sec_kind == "Round tube":
    d_ = sb.number_input("Outer diameter d (mm)", 10.0, 1000.0, 100.0)
    t_ = sb.number_input("Wall t (mm)", 1.0, 200.0, 5.0)
    sec = bc.tube(d_, t_)
    if sec is None:
        st.error("Wall thickness is too large for this diameter.")
        st.stop()
else:
    sec = bc.SectionProps(sb.number_input("I (mm4)", 1e3, 1e11, 8e7, format="%.4g"),
                          sb.number_input("c (mm)", 1.0, 2000.0, 100.0), None)
if mass_per_m is None and sec.area_mm2:
    mass_per_m = sec.area_mm2 * 1e-6 * density
use_sw = sb.checkbox("Include self-weight", value=True)
sw = bc.self_weight(mass_per_m) if (use_sw and mass_per_m) else 0.0
limit_den = sb.select_slider("Deflection limit L /", options=[180, 250, 300, 360, 400, 500], value=360)
safety = sb.number_input("Safety factor", 1.0, 5.0, 1.5, 0.1)

# ----------------------------------------------------------------- main: geometry and loads
left, right = st.columns([1, 1])
with left:
    st.subheader("Geometry")
    if mode == "Single span":
        supports = st.selectbox("Supports", list(bc.SUPPORTS), format_func=lambda s: {
            "simple": "Simply supported", "cantilever": "Cantilever (fixed at x = 0)",
            "fixed": "Fixed both ends", "propped": "Propped cantilever (fixed at x = 0)"}[s])
        length = st.number_input("Span L (m)", 0.5, 50.0, 5.0, 0.5)
    else:
        spans_txt = st.text_input("Span lengths (m), comma separated", "4, 4, 4")
        try:
            spans = [float(s) for s in spans_txt.replace(";", ",").split(",") if s.strip()]
            assert spans and all(s > 0 for s in spans)
        except Exception:
            st.error("Enter positive numbers such as 4, 5, 4")
            st.stop()
        overhang_l = st.number_input("Left overhang (m)", 0.0, 20.0, 0.0, 0.5)
        overhang_r = st.number_input("Right overhang (m)", 0.0, 20.0, 0.0, 0.5)
        clamp_l = st.checkbox("Clamp left end (fixed)", value=False)
        clamp_r = st.checkbox("Clamp right end (fixed)", value=False)
        sup_x = [overhang_l]
        for s in spans:
            sup_x.append(sup_x[-1] + s)
        length = sup_x[-1] + overhang_r
        st.write(f"Total length {length:g} m, supports at {', '.join(f'{x:g}' for x in sup_x)} m")

with right:
    st.subheader("Loads")
    udl_kn = st.number_input("Uniform load over the whole length (kN/m)", 0.0, 1000.0, 6.0, 0.5)
    pl_df = st.data_editor(pd.DataFrame({"load_kN": [20.0], "pos_m": [min(2.0, length)]}), num_rows="dynamic",
                           key="pl", use_container_width=True)
    pu_df = st.data_editor(pd.DataFrame({"w_kN_per_m": [0.0], "start_m": [0.0], "end_m": [min(2.0, length)]}),
                           num_rows="dynamic", key="pu", use_container_width=True)

point_loads = [bc.PointLoad(float(r.load_kN) * 1e3, float(r.pos_m)) for r in pl_df.dropna().itertuples()
               if r.load_kN > 0]
partials = [bc.PartialUdl(float(r.w_kN_per_m) * 1e3, float(r.start_m), float(r.end_m)) for r in pu_df.dropna().itertuples()
            if r.w_kN_per_m > 0 and r.end_m > r.start_m]


# ----------------------------------------------------------------- analysis
def analyse(e):
    """Return dict of arrays and numbers for modulus e (GPa)."""
    if mode == "Single span":
        g = bc.GeneralInput(supports, length, e, sec.i_mm4, point_loads, udl_kn * 1e3, partials, sw, sec.c_mm,
                            limit_den, strength, safety, sec.shear_area_mm2, sec.shear_factor, shear_strength)
        r = bc.solve(g)
        p = bc.curves(g, 301)
        return {"x": p.x, "d": p.deflection_mm, "m": p.moment_nm, "v": p.shear_n, "reactions": None,
                "dmax": r.deflection_mm, "dat": r.deflection_at_m, "sag": r.sag_moment_nm, "hog": r.hog_moment_nm,
                "vmax": r.shear_n, "stress": r.stress_mpa, "allow": r.allowable_mpa, "util": r.utilisation,
                "ratio": r.span_ratio, "ok": r.deflection_ok, "shear_util": r.shear_utilisation, "inp": g, "res": r}
    c = bc.ContinuousInput(length, sup_x, e, sec.i_mm4, point_loads, udl_kn * 1e3, partials, sw, clamp_l, clamp_r,
                           sec.c_mm, strength, safety, limit_den)
    r = bc.solve_continuous(c)
    return {"x": r.x, "d": r.deflection_mm, "m": r.moment_nm, "v": r.shear_n, "reactions": r.reactions,
            "dmax": r.deflection_max_mm, "dat": r.deflection_at_m, "sag": r.sag_moment_nm, "hog": r.hog_moment_nm,
            "vmax": r.shear_max_n, "stress": r.stress_mpa, "allow": r.allowable_mpa, "util": r.utilisation,
            "ratio": r.span_ratio, "ok": r.deflection_ok, "shear_util": None, "inp": c, "res": r}


try:
    out = analyse(e_gpa)
except ValueError as err:
    st.error(str(err))
    st.stop()

tab_res, tab_diag, tab_3d, tab_pat, tab_tools = st.tabs(["Results", "Diagrams", "3D view", "Load patterns", "Tools"])

with tab_res:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Max deflection", f"{out['dmax']:.2f} mm", f"at x = {out['dat']:.2f} m", delta_color="off")
    c2.metric("Span / deflection", "inf" if math.isinf(out["ratio"]) else f"{out['ratio']:.0f}",
              "OK" if out["ok"] else "exceeds limit", delta_color="normal" if out["ok"] else "inverse")
    c3.metric("Max moment", f"{max(out['sag'], -out['hog']) / 1e3:.2f} kNm")
    c4.metric("Max shear", f"{out['vmax'] / 1e3:.2f} kN")
    if out["util"] is not None:
        u1, u2 = st.columns(2)
        u1.metric("Bending stress", f"{out['stress']:.1f} MPa", f"allowable {out['allow']:.1f} MPa", delta_color="off")
        u2.metric("Bending utilisation", f"{out['util']:.2f}", "OK" if out["util"] <= 1 else "FAIL",
                  delta_color="normal" if out["util"] <= 1 else "inverse")
    if out["shear_util"] is not None:
        st.write(f"Shear utilisation: {out['shear_util']:.2f}")
    if out["reactions"]:
        rows = [{"position (m)": (k if not isinstance(k, tuple) else k[0]),
                 "type": "reaction (kN)" if not isinstance(k, tuple) else "clamp moment (kNm)",
                 "value": v / 1e3} for k, v in out["reactions"].items()]
        st.dataframe(pd.DataFrame(rows).round(3), use_container_width=True)
        st.pyplot(bp.plot_reactions(out["reactions"], length))

with tab_diag:
    sx = sup_x if mode == "Continuous beam" else None
    st.pyplot(bp.plot_profiles(out["x"], out["d"], out["m"], out["v"], sx))

with tab_3d:
    st.pyplot(bp.plot_beam_3d(out["x"], out["d"], out["m"], length))

with tab_pat:
    if mode != "Continuous beam":
        st.info("Live-load patterns apply to continuous beams. Switch the beam type in the sidebar.")
    else:
        st.write("The permanent loads above (uniform load, tables, self-weight) are factored with gamma_G. "
                 "The live load below is applied on every combination of spans (2^n patterns).")
        p1, p2, p3, p4 = st.columns(4)
        dead_extra = p1.number_input("Extra permanent (kN/m)", 0.0, 1000.0, 0.0, 0.5)
        live_kn = p2.number_input("Live load Q (kN/m)", 0.0, 1000.0, 10.0, 0.5)
        g_g = p3.number_input("gamma G", 1.0, 2.0, 1.35, 0.05)
        g_q = p4.number_input("gamma Q", 1.0, 2.5, 1.5, 0.05)
        try:
            uls_env = bc.pattern_envelope(out["inp"], dead_extra * 1e3, live_kn * 1e3, g_g, g_q)
            sls_env = bc.pattern_envelope(out["inp"], dead_extra * 1e3, live_kn * 1e3, 1.0, 1.0)
        except ValueError as err:
            st.error(str(err))
        else:
            m_uls = max(uls_env["m_max"].max(), -uls_env["m_min"].min())
            q1, q2, q3, q4 = st.columns(4)
            q1.metric("Patterns tried", len(uls_env["patterns"]))
            q2.metric("ULS max moment", f"{m_uls / 1e3:.1f} kNm")
            q3.metric("ULS max shear", f"{max(uls_env['v_max'].max(), -uls_env['v_min'].min()) / 1e3:.1f} kN")
            q4.metric("SLS max deflection", f"{sls_env['d_max'].max():.1f} mm")
            if strength:
                util_u = m_uls * 1000.0 * sec.c_mm / sec.i_mm4 / (strength / safety)
                st.write(f"ULS bending utilisation (strength / safety factor): **{util_u:.2f}**")
            st.pyplot(bp.plot_envelopes(uls_env, sup_x, "ULS envelopes"))
            rk = pd.DataFrame({"Support x (m)": uls_env["r_keys"], "Min reaction (kN)": uls_env["r_min"] / 1e3,
                               "Max reaction (kN)": uls_env["r_max"] / 1e3}).round(2)
            st.dataframe(rk, use_container_width=True)

with tab_tools:
    st.subheader("Material comparison (same section and loads)")
    rows = []
    for m_ in bc.MATERIALS:
        o = analyse(m_.e_gpa)
        rows.append({"Material": m_.name, "E (GPa)": m_.e_gpa, "Deflection (mm)": round(o["dmax"], 2),
                     "L/defl": "inf" if math.isinf(o["ratio"]) else round(o["ratio"])})
    st.dataframe(pd.DataFrame(rows), use_container_width=True)

    if mode == "Single span":
        st.subheader("Lightest passing IPE")
        sel = bc.select_ipe(supports, length, bc.Material("sel", e_gpa, density, strength, shear_strength),
                            point_loads, udl_kn * 1e3 + 0.0, limit_den, safety, use_sw, n=151)
        sel_df = pd.DataFrame(sel).round(2)
        st.dataframe(sel_df, use_container_width=True)
        first = next((r["section"] for r in sel if r["ok"]), None)
        st.write(f"Lightest passing: **{first}**" if first else "No IPE size passes all checks.")
        if mass_per_m:
            f1 = bc.natural_frequency(supports, length, e_gpa, sec.i_mm4, mass_per_m)
            st.write(f"First natural frequency (beam mass only): **{f1:.1f} Hz**")

    st.subheader("Download")
    prof = pd.DataFrame({"x_m": out["x"], "deflection_mm": out["d"], "moment_Nm": out["m"], "shear_N": out["v"]})
    st.download_button("Profile CSV", prof.round(4).to_csv(index=False), "beam_profile.csv", "text/csv")
    if mode == "Single span":
        st.download_button("Report (markdown)", bc.report_text(out["inp"], out["res"], "Beam report"),
                           "beam_report.md", "text/markdown")

    # PDF report
    lines = [f"Beam type : {mode}", f"Length    : {length:g} m", f"Material  : {mat_name} (E = {e_gpa} GPa)",
             f"Section   : {sec_kind}, I = {sec.i_mm4:.4g} mm4", f"Self-weight: {sw:.1f} N/m",
             f"Uniform load: {udl_kn:g} kN/m", "", f"Max deflection : {out['dmax']:.2f} mm at x = {out['dat']:.2f} m",
             f"Span / deflection: {'inf' if math.isinf(out['ratio']) else format(out['ratio'], '.0f')} (limit L/{limit_den})",
             f"Max sagging moment : {out['sag'] / 1e3:.2f} kNm", f"Max hogging moment : {out['hog'] / 1e3:.2f} kNm",
             f"Max shear : {out['vmax'] / 1e3:.2f} kN"]
    if out["util"] is not None:
        lines.append(f"Bending utilisation : {out['util']:.2f}")
    buf = io.BytesIO()
    bp.make_pdf_report(buf, "BeamDeflection report", lines, out["x"], out["d"], out["m"], out["v"], length,
                       sup_x if mode == "Continuous beam" else None)
    st.download_button("PDF report", buf.getvalue(), "beam_report.pdf", "application/pdf")
