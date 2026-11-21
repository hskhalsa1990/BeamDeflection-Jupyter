"""Run: python test_beam_core.py  (or pytest)"""
import math
from beam_core import *

E, I, L = 210.0, 8.0e7, 4.0
EI = E * 1e9 * I * 1e-12


def close(a, b, tol=1e-3):
    assert abs(a - b) <= tol * max(1.0, abs(b)), (a, b)


def test_simple_point_centre():
    r = solve(GeneralInput("simple", L, E, I, [PointLoad(10000, L / 2)]))
    close(r.deflection_mm, 10000 * L ** 3 / (48 * EI) * 1000)
    close(r.max_moment_nm, 10000 * L / 4)


def test_simple_udl():
    r = solve(GeneralInput("simple", L, E, I, udl_n_per_m=5000))
    close(r.deflection_mm, 5 * 5000 * L ** 4 / (384 * EI) * 1000)


def test_cantilever_tip():
    r = solve(GeneralInput("cantilever", L, E, I, [PointLoad(2000, L)]))
    close(r.deflection_mm, 2000 * L ** 3 / (3 * EI) * 1000)
    close(r.hog_moment_nm, -2000 * L)


def test_fixed_udl_and_partial_equals_full():
    full = solve(GeneralInput("fixed", L, E, I, udl_n_per_m=3000))
    part = solve(GeneralInput("fixed", L, E, I, partial_udls=[PartialUdl(3000, 0, L)]))
    close(full.deflection_mm, 3000 * L ** 4 / (384 * EI) * 1000)
    close(part.deflection_mm, full.deflection_mm)


def test_propped_udl():
    r = solve(GeneralInput("propped", L, E, I, udl_n_per_m=4000))
    close(r.deflection_mm, 0.0054161216 * 4000 * L ** 4 / EI * 1000, 2e-3)


def test_general_matches_standard_cases():
    for sup in SUPPORTS:
        for uni in (False, True):
            c = standard_case(sup, uni)
            if uni:
                g = solve(GeneralInput(sup, L, E, I, udl_n_per_m=1000))
                d, m, v = quick_max(c, 1000, L, E, I)
            else:
                pos = L if sup == "cantilever" else L / 2
                g = solve(GeneralInput(sup, L, E, I, [PointLoad(1000, pos)]))
                d, m, v = quick_max(c, 1000, L, E, I)
            close(g.deflection_mm, d, 3e-3)
            close(g.max_moment_nm, m, 3e-3)
            close(g.shear_n, v, 3e-3)


def test_select_ipe_monotonic():
    rows = select_ipe("simple", 5.0, material("S275"), udl_n_per_m=10000)
    defl = [r["deflection_mm"] for r in rows]
    assert all(a > b for a, b in zip(defl, defl[1:]))
    assert lightest_ipe("simple", 5.0, material("S275"), udl_n_per_m=10000) is not None


def test_influence_symmetry():
    xs = [0.25, 0.75]
    d = influence_deflection("simple", xs, [0.25, 0.75])
    close(d[0][1], d[1][0])   # Maxwell reciprocity


def test_natural_frequency_simple():
    # f1 = (pi / 2L^2) sqrt(EI/m)
    f = natural_frequency("simple", L, E, I, 50.0)
    close(f, math.pi / (2 * L ** 2) * math.sqrt(EI / 50.0))
    assert natural_frequency("fixed", L, E, I, 50.0) > natural_frequency("simple", L, E, I, 50.0) > \
        natural_frequency("cantilever", L, E, I, 50.0)


def test_moving_load_envelope():
    env = moving_load_envelope("simple", L, 10000, E, I, n_pos=41, n_x=81)
    close(env["sag_env"].max(), 10000 * L / 4)          # worst moment: load at midspan
    close(env["shear_env"].max(), 10000, 1e-6)         # load at the support


def test_continuous_two_span_udl():
    w = 5000.0
    r = solve_continuous(ContinuousInput(2 * L, equal_spans(2, L), E, I, udl_n_per_m=w))
    close(r.reactions[L], 1.25 * w * L)
    close(r.reactions[0.0], 0.375 * w * L)
    close(r.hog_moment_nm, -w * L * L / 8)
    close(r.sag_moment_nm, 9 * w * L * L / 128)
    close(r.deflection_max_mm, 0.0054161216 * w * L ** 4 / EI * 1000, 2e-3)


def test_continuous_matches_single_span_and_cantilever():
    r = solve_continuous(ContinuousInput(L, [0.0, L], E, I, point_loads=[PointLoad(1e4, L / 2)]))
    close(r.deflection_max_mm, 1e4 * L ** 3 / (48 * EI) * 1000)
    c = solve_continuous(ContinuousInput(L, [], E, I, point_loads=[PointLoad(2e3, L)], fixed_left=True))
    close(c.deflection_max_mm, 2e3 * L ** 3 / (3 * EI) * 1000)
    close(c.hog_moment_nm, -2e3 * L)


def test_continuous_unstable_raises():
    try:
        solve_continuous(ContinuousInput(L, [0.0], E, I, udl_n_per_m=1000))
    except ValueError:
        return
    raise AssertionError("expected ValueError")


def test_pattern_envelope_two_span():
    w = 5000.0
    env = pattern_envelope(ContinuousInput(2 * L, equal_spans(2, L), E, I), 0.0, w, 1.0, 1.0)
    close(env["m_min"].min(), -w * L * L / 8)                 # both spans loaded
    close(env["m_max"].max(), 49.0 / 512.0 * w * L * L, 2e-3)  # one span loaded
    assert len(env["patterns"]) == 4


def test_load_combination():
    close(load_combination(10.0, 5.0), 1.35 * 10 + 1.5 * 5)


def test_shear_deflection():
    g = 210e9 / (2.6)
    close(shear_deflection_mm("simple", True, 1000.0, 4.0, 210.0, 3000.0), 1000 * 16 / (8 * g * 3e-3) * 1000)
    close(shear_deflection_mm("cantilever", False, 1000.0, 4.0, 210.0, 3000.0), 1000 * 4 / (g * 3e-3) * 1000)
    assert shear_deflection_mm("propped", True, 1.0, 1.0, 210.0, 1.0) is None


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            f()
            print("ok", k)
