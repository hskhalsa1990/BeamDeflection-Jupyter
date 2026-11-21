"""beam_plots.py - matplotlib figures shared by the notebook helpers and the Streamlit dashboard."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def plot_profiles(x, defl_mm, moment_nm, shear_n, supports_x=None, title=""):
    """Three stacked diagrams: deflection (downward negative), moment (sagging positive), shear."""
    fig, ax = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True)
    spec = [(-np.asarray(defl_mm), "Deflection (mm)", "#1f77b4"),
            (np.asarray(moment_nm) / 1e3, "Moment (kNm)", "#d62728"),
            (np.asarray(shear_n) / 1e3, "Shear (kN)", "#2ca02c")]
    for a, (y, lab, col) in zip(ax, spec):
        a.plot(x, y, color=col, lw=2)
        a.fill_between(x, y, 0, color=col, alpha=0.15)
        a.axhline(0, color="k", lw=0.8)
        a.set_ylabel(lab)
        a.grid(alpha=0.25)
        for sx in (supports_x or []):
            a.axvline(sx, color="grey", lw=0.8, ls=":")
    ax[-1].set_xlabel("x (m)")
    if title:
        ax[0].set_title(title)
    fig.tight_layout()
    return fig


def plot_beam_3d(x, defl_mm, moment_nm, length_m, width_m=None, title="3D deformed beam (colour = moment)"):
    """Deformed beam surface coloured by bending moment; deflection exaggerated for visibility."""
    x = np.asarray(x); d = np.asarray(defl_mm); m = np.asarray(moment_nm) / 1e3
    width_m = width_m or 0.08 * length_m
    ys = np.linspace(-width_m / 2, width_m / 2, 7)
    X, Y = np.meshgrid(x, ys)
    k = 0.12 * length_m / max(float(np.abs(d).max()), 1e-9)
    Z = np.tile(-d * k, (len(ys), 1))
    M = np.tile(m, (len(ys), 1))
    norm = Normalize(M.min(), M.max() if M.max() > M.min() else M.min() + 1)
    fig = plt.figure(figsize=(10, 5.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, facecolors=cm.coolwarm(norm(M)), rstride=1, cstride=1, antialiased=True)
    ax.plot_surface(X, Y, np.zeros_like(Z), color="grey", alpha=0.1)
    sm = cm.ScalarMappable(norm=norm, cmap="coolwarm"); sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.08, label="Bending moment (kNm)")
    ax.set_xlabel("x (m)"); ax.set_ylabel("width (m)"); ax.set_zlabel("deflection (exaggerated)")
    ax.set_title(title); ax.view_init(elev=24, azim=-62)
    ax.set_box_aspect((3, 1, 0.8))
    fig.tight_layout()
    return fig


def plot_reactions(reactions, length_m):
    """Bar chart of support reactions (kN)."""
    items = [(k, v) for k, v in reactions.items() if not isinstance(k, tuple)]
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.bar([f"{k:g} m" for k, _ in items], [v / 1e3 for _, v in items], color="#1f77b4")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("Reaction (kN, up +)")
    ax.set_title("Support reactions")
    fig.tight_layout()
    return fig


def plot_envelopes(env, supports_x=None, title="Load-pattern envelopes"):
    """Min/max envelopes of moment, shear and deflection from beam_core.pattern_envelope."""
    x = env["x"]
    fig, ax = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    spec = [(env["m_max"] / 1e3, env["m_min"] / 1e3, "Moment (kNm)", "#d62728"),
            (env["v_max"] / 1e3, env["v_min"] / 1e3, "Shear (kN)", "#2ca02c"),
            (-env["d_min"], -env["d_max"], "Deflection (mm)", "#1f77b4")]
    for a, (hi, lo, lab, col) in zip(ax, spec):
        a.fill_between(x, lo, hi, color=col, alpha=0.2)
        a.plot(x, hi, color=col, lw=1.6); a.plot(x, lo, color=col, lw=1.6)
        a.axhline(0, color="k", lw=0.8); a.set_ylabel(lab); a.grid(alpha=0.25)
        for sx in (supports_x or []):
            a.axvline(sx, color="grey", lw=0.8, ls=":")
    ax[-1].set_xlabel("x (m)"); ax[0].set_title(title)
    fig.tight_layout()
    return fig


def make_pdf_report(target, title, lines, x, defl_mm, moment_nm, shear_n, length_m, supports_x=None):
    """A4 PDF: page 1 inputs and results, page 2 diagrams, page 3 3D beam. target = path or file-like object."""
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(target) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.08, 0.95, title, fontsize=16, weight="bold")
        y = 0.91
        for ln in lines:
            fig.text(0.08, y, ln, fontsize=9.5, family="monospace", va="top")
            y -= 0.022
            if y < 0.05:
                break
        pdf.savefig(fig); plt.close(fig)
        for f in (plot_profiles(x, defl_mm, moment_nm, shear_n, supports_x, "Deflection, moment, shear"),
                  plot_beam_3d(x, defl_mm, moment_nm, length_m)):
            pdf.savefig(f); plt.close(f)
