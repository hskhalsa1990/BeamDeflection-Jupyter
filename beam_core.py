"""beam_core.py - Euler-Bernoulli beam deflection, moment and shear.

Python port of the BeamDeflection Android app (Kotlin domain + data layers).
Single span, constant section, linear-elastic, small deflections, strong-axis bending.
Units: lengths m (sections in mm), forces N, loads N/m, E in GPa, I in mm^4.
"""
from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

GRAVITY = 9.80665
SUPPORTS = ("simple", "cantilever", "fixed", "propped")


# --------------------------------------------------------------------------
# Data: materials and sections
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Material:
    name: str
    e_gpa: float
    density: float  # kg/m3
    strength_mpa: Optional[float]
    shear_strength_mpa: Optional[float] = None


MATERIALS = [
    Material("Structural steel S235", 210.0, 7850.0, 235.0, 141.0),
    Material("Structural steel S275", 210.0, 7850.0, 275.0, 165.0),
    Material("Structural steel S355", 210.0, 7850.0, 355.0, 213.0),
    Material("Stainless steel 304", 193.0, 8000.0, 215.0, 129.0),
    Material("Aluminium 6061-T6", 68.9, 2700.0, 276.0, 166.0),
    Material("Aluminium 6063-T5", 68.9, 2700.0, 145.0, 87.0),
    Material("Titanium Ti-6Al-4V", 113.8, 4430.0, 880.0, 528.0),
    Material("Copper (annealed)", 117.0, 8960.0, 70.0, 42.0),
    Material("Grey cast iron", 110.0, 7200.0, None),
    Material("Concrete C30/37", 33.0, 2400.0, None),
    Material("Timber C24 (softwood)", 11.0, 420.0, 24.0, 4.0),
    Material("Timber D40 (hardwood)", 11.0, 660.0, 40.0, 4.0),
    Material("Glulam GL24h", 11.5, 420.0, 24.0, 3.5),
    Material("CFRP laminate (quasi-isotropic)", 70.0, 1600.0, None),
]


def material(name_part: str) -> Material:
    for m in MATERIALS:
        if name_part.lower() in m.name.lower():
            return m
    raise KeyError(name_part)


@dataclass(frozen=True)
class SectionProps:
    i_mm4: float
    c_mm: float
    area_mm2: Optional[float]
    mass_kg_per_m: Optional[float] = None
    shear_area_mm2: Optional[float] = None
    shear_factor: float = 1.0


# name, h mm, Iy cm4, Wy cm3, mass kg/m, web thickness mm
IPE = [
    ("IPE 100", 100.0, 171.0, 34.2, 8.1, 4.1),
    ("IPE 120", 120.0, 318.0, 53.0, 10.4, 4.4),
    ("IPE 140", 140.0, 541.0, 77.3, 12.9, 4.7),
    ("IPE 160", 160.0, 869.0, 109.0, 15.8, 5.0),
    ("IPE 180", 180.0, 1317.0, 146.0, 18.8, 5.3),
    ("IPE 200", 200.0, 1943.0, 194.0, 22.4, 5.6),
    ("IPE 220", 220.0, 2772.0, 252.0, 26.2, 5.9),
    ("IPE 240", 240.0, 3892.0, 324.0, 30.7, 6.2),
    ("IPE 270", 270.0, 5790.0, 429.0, 36.1, 6.6),
    ("IPE 300", 300.0, 8356.0, 557.0, 42.2, 7.1),
]


def rect(b: float, h: float) -> SectionProps:
    return SectionProps(b * h ** 3 / 12.0, h / 2.0, b * h, shear_area_mm2=b * h, shear_factor=1.5)


def round_bar(d: float) -> SectionProps:
    a = math.pi * d * d / 4.0
    return SectionProps(math.pi * d ** 4 / 64.0, d / 2.0, a, shear_area_mm2=a, shear_factor=4.0 / 3.0)


def box(bo: float, ho: float, t: float) -> Optional[SectionProps]:
    bi, hi = bo - 2 * t, ho - 2 * t
    if bi <= 0 or hi <= 0:
        return None
    return SectionProps((bo * ho ** 3 - bi * hi ** 3) / 12.0, ho / 2.0, bo * ho - bi * hi,
                        shear_area_mm2=2.0 * t * ho)


def tube(d: float, t: float) -> Optional[SectionProps]:
    di = d - 2 * t
    if di <= 0:
        return None
    a = math.pi * (d * d - di * di) / 4.0
    return SectionProps(math.pi * (d ** 4 - di ** 4) / 64.0, d / 2.0, a, shear_area_mm2=a,
                        shear_factor=4.0 / 3.0 * (d * d + d * di + di * di) / (d * d + di * di))


def ipe(name: str) -> SectionProps:
    for n, h, iy, _wy, mass, tw in IPE:
        if n == name:
            return SectionProps(iy * 1e4, h / 2.0, None, mass, shear_area_mm2=h * tw)
    raise KeyError(name)


def self_weight(mass_kg_per_m: Optional[float]) -> float:
    return (mass_kg_per_m or 0.0) * GRAVITY


def self_weight_from_area(area_mm2: float, density: float) -> float:
    return area_mm2 * 1e-6 * density * GRAVITY


# --------------------------------------------------------------------------
# Standard single-load cases (closed form, vectorised over xi = x / L)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class LoadCase:
    key: str
    label: str
    uniform: bool
    supports: str
    defl_coeff: float
    moment_coeff: float
    shear_coeff: float
    formula: str

    @property
    def span_power(self) -> int:
        return 4 if self.uniform else 3

    @property
    def partner(self) -> "LoadCase":
        return CASES[(self.supports, True)]

    def shape(self, xi):
        xi = np.asarray(xi, dtype=float)
        s = np.minimum(xi, 1.0 - xi)
        k = (self.supports, self.uniform)
        if k == ("simple", False):
            return s * (3.0 - 4.0 * s * s) / 48.0
        if k == ("simple", True):
            return xi * (1.0 - 2.0 * xi ** 2 + xi ** 3) / 24.0
        if k == ("cantilever", False):
            return xi ** 2 * (3.0 - xi) / 6.0
        if k == ("cantilever", True):
            return xi ** 2 * (6.0 - 4.0 * xi + xi ** 2) / 24.0
        if k == ("fixed", False):
            return s * s * (3.0 - 4.0 * s) / 48.0
        if k == ("fixed", True):
            return xi ** 2 * (1.0 - xi) ** 2 / 24.0
        if k == ("propped", False):
            return np.where(xi <= 0.5, xi ** 2 * (9.0 - 11.0 * xi) / 96.0,
                            (-2.0 + 12.0 * xi - 15.0 * xi ** 2 + 5.0 * xi ** 3) / 96.0)
        return xi ** 2 * (1.0 - xi) * (3.0 - 2.0 * xi) / 48.0

    def moment_shape(self, xi):
        xi = np.asarray(xi, dtype=float)
        s = np.minimum(xi, 1.0 - xi)
        k = (self.supports, self.uniform)
        if k == ("simple", False):
            return s / 2.0
        if k == ("simple", True):
            return xi * (1.0 - xi) / 2.0
        if k == ("cantilever", False):
            return -(1.0 - xi)
        if k == ("cantilever", True):
            return -(1.0 - xi) ** 2 / 2.0
        if k == ("fixed", False):
            return s / 2.0 - 1.0 / 8.0
        if k == ("fixed", True):
            return xi * (1.0 - xi) / 2.0 - 1.0 / 12.0
        if k == ("propped", False):
            return np.where(xi <= 0.5, 11.0 * xi / 16.0 - 3.0 / 16.0, 5.0 * (1.0 - xi) / 16.0)
        return 5.0 * xi / 8.0 - xi ** 2 / 2.0 - 1.0 / 8.0

    def shear_shape(self, xi):
        xi = np.asarray(xi, dtype=float)
        k = (self.supports, self.uniform)
        if k in (("simple", False), ("fixed", False)):
            return np.where(xi < 0.5, 0.5, -0.5)
        if k in (("simple", True), ("fixed", True)):
            return (1.0 - 2.0 * xi) / 2.0
        if k == ("cantilever", False):
            return np.ones_like(xi)
        if k == ("cantilever", True):
            return 1.0 - xi
        if k == ("propped", False):
            return np.where(xi < 0.5, 11.0 / 16.0, -5.0 / 16.0)
        return 5.0 / 8.0 - xi


def _mk():
    d = {}
    rows = [
        ("simple", False, "Simply supported - point load at centre", 1 / 48, 1 / 4, 1 / 2, "d = PL^3 / 48EI"),
        ("simple", True, "Simply supported - uniform load", 5 / 384, 1 / 8, 1 / 2, "d = 5wL^4 / 384EI"),
        ("cantilever", False, "Cantilever - point load at free end", 1 / 3, 1.0, 1.0, "d = PL^3 / 3EI"),
        ("cantilever", True, "Cantilever - uniform load", 1 / 8, 1 / 2, 1.0, "d = wL^4 / 8EI"),
        ("fixed", False, "Fixed both ends - point load at centre", 1 / 192, 1 / 8, 1 / 2, "d = PL^3 / 192EI"),
        ("fixed", True, "Fixed both ends - uniform load", 1 / 384, 1 / 12, 1 / 2, "d = wL^4 / 384EI"),
        ("propped", False, "Propped cantilever - point load at centre", 1 / (48 * math.sqrt(5)), 3 / 16, 11 / 16,
         "d ~ 0.00932 PL^3 / EI"),
        ("propped", True, "Propped cantilever - uniform load", 0.0054161216, 1 / 8, 5 / 8, "d ~ 0.00542 wL^4 / EI"),
    ]
    for sup, uni, label, dc, mc, sc, f in rows:
        d[(sup, uni)] = LoadCase(f"{sup}_{'udl' if uni else 'point'}", label, uni, sup, dc, mc, sc, f)
    return d


CASES = _mk()


def standard_case(supports: str, uniform: bool) -> LoadCase:
    return CASES[(supports, uniform)]


def quick_max(case: LoadCase, load: float, length_m: float, e_gpa: float, i_mm4: float):
    """Closed-form max deflection (mm), max moment (N m) and max shear/reaction (N)."""
    ei = e_gpa * 1e9 * i_mm4 * 1e-12
    defl = case.defl_coeff * load * length_m ** case.span_power / ei * 1000.0
    if case.uniform:
        return defl, case.moment_coeff * load * length_m ** 2, case.shear_coeff * load * length_m
    return defl, case.moment_coeff * load * length_m, case.shear_coeff * load


# --------------------------------------------------------------------------
# General beam: point loads anywhere + uniform loads (full / partial) + self weight
# --------------------------------------------------------------------------
@dataclass
class PointLoad:
    load_n: float
    pos_m: float


@dataclass
class PartialUdl:
    w_n_per_m: float
    start_m: float
    end_m: float


@dataclass
class GeneralInput:
    supports: str
    length_m: float
    e_gpa: float
    i_mm4: float
    point_loads: List[PointLoad] = field(default_factory=list)
    udl_n_per_m: float = 0.0
    partial_udls: List[PartialUdl] = field(default_factory=list)
    self_weight_n_per_m: float = 0.0
    c_mm: Optional[float] = None
    limit_denominator: int = 360
    strength_mpa: Optional[float] = None
    safety_factor: float = 1.5
    shear_area_mm2: Optional[float] = None
    shear_factor: float = 1.0
    shear_strength_mpa: Optional[float] = None


@dataclass
class GeneralResult:
    deflection_mm: float
    deflection_at_m: float
    sag_moment_nm: float
    hog_moment_nm: float
    max_moment_nm: float
    shear_n: float
    stress_mpa: Optional[float]
    allowable_mpa: Optional[float]
    utilisation: Optional[float]
    shear_stress_mpa: Optional[float]
    shear_allowable_mpa: Optional[float]
    shear_utilisation: Optional[float]
    span_ratio: float
    deflection_ok: bool


@dataclass
class Profile:
    x: np.ndarray
    deflection_mm: np.ndarray
    moment_nm: np.ndarray
    shear_n: np.ndarray


def _propped_reaction(p, l, a):
    return p * a * a * (3.0 * l - a) / (2.0 * l ** 3)


def _fixed_reaction_a(p, l, a):
    b = l - a
    return p * b * b * (3.0 * a + b) / l ** 3


def _cantilever_y(p, a, x):
    if x <= a:
        return p * x * x * (3.0 * a - x) / 6.0
    return p * a * a * (3.0 * x - a) / 6.0


def _point_ei_y(s, p, l, a, x):
    b = l - a
    if s == "simple":
        if x <= a:
            return p * b * x * (l * l - b * b - x * x) / (6.0 * l)
        xr = l - x
        return p * a * xr * (l * l - a * a - xr * xr) / (6.0 * l)
    if s == "cantilever":
        return _cantilever_y(p, a, x)
    if s == "fixed":
        if x <= a:
            return p * b * b * x * x * (3.0 * a * l - x * (3.0 * a + b)) / (6.0 * l ** 3)
        xr = l - x
        return p * a * a * xr * xr * (3.0 * b * l - xr * (3.0 * b + a)) / (6.0 * l ** 3)
    return _cantilever_y(p, a, x) - _propped_reaction(p, l, a) * x * x * (3.0 * l - x) / 6.0


def _point_moment(s, p, l, a, x):
    b = l - a
    if s == "simple":
        return p * b * x / l if x <= a else p * a * (l - x) / l
    if s == "cantilever":
        return -p * (a - x) if x < a else 0.0
    if s == "fixed":
        ra = _fixed_reaction_a(p, l, a)
        ma = -p * a * b * b / (l * l)
        return ra * x + ma - (p * (x - a) if x > a else 0.0)
    return _propped_reaction(p, l, a) * (l - x) - (p * (a - x) if x < a else 0.0)


def _point_shear(s, p, l, a, x):
    b = l - a
    before = x < a
    if s == "simple":
        ra = p * b / l
        return ra if before else ra - p
    if s == "cantilever":
        return p if before else 0.0
    if s == "fixed":
        ra = _fixed_reaction_a(p, l, a)
        return ra if before else ra - p
    rb = _propped_reaction(p, l, a)
    return p - rb if before else -rb


_GL_T = 0.7745966692414834


def _gauss3(lo, hi, f):
    if hi <= lo:
        return 0.0
    mid, half = (lo + hi) / 2.0, (hi - lo) / 2.0
    return half * (5.0 / 9.0 * f(mid - half * _GL_T) + 8.0 / 9.0 * f(mid) + 5.0 / 9.0 * f(mid + half * _GL_T))


def _over_load(lo, hi, x, f):
    return _gauss3(lo, min(hi, x), f) + _gauss3(max(lo, x), hi, f)


def _clean_loads(inp: GeneralInput):
    l = inp.length_m
    pts = [PointLoad(p.load_n, min(max(p.pos_m, 0.0), l)) for p in inp.point_loads if p.load_n > 0.0]
    udls = [PartialUdl(u.w_n_per_m, min(max(u.start_m, 0.0), l), min(max(u.end_m, 0.0), l))
            for u in inp.partial_udls]
    udls = [u for u in udls if u.w_n_per_m > 0.0 and u.end_m > u.start_m]
    return pts, udls


def profile(inp: GeneralInput, xs) -> Profile:
    l = inp.length_m
    ei = inp.e_gpa * 1e9 * inp.i_mm4 * 1e-12
    pts, udls = _clean_loads(inp)
    w = inp.udl_n_per_m + inp.self_weight_n_per_m
    uc = standard_case(inp.supports, True)
    xs = np.asarray(xs, dtype=float)
    d = np.zeros(len(xs))
    m = np.zeros(len(xs))
    v = np.zeros(len(xs))
    s = inp.supports
    for i, x in enumerate(xs):
        y = mm = vv = 0.0
        for pl in pts:
            y += _point_ei_y(s, pl.load_n, l, pl.pos_m, x) / ei
            mm += _point_moment(s, pl.load_n, l, pl.pos_m, x)
            vv += _point_shear(s, pl.load_n, l, pl.pos_m, x)
        for u in udls:
            y += u.w_n_per_m * _over_load(u.start_m, u.end_m, x, lambda a: _point_ei_y(s, 1.0, l, a, x)) / ei
            mm += u.w_n_per_m * _over_load(u.start_m, u.end_m, x, lambda a: _point_moment(s, 1.0, l, a, x))
            vv += u.w_n_per_m * _over_load(u.start_m, u.end_m, x, lambda a: _point_shear(s, 1.0, l, a, x))
        if w > 0.0:
            xi = x / l
            y += w * l ** 4 / ei * float(uc.shape(xi))
            mm += w * l * l * float(uc.moment_shape(xi))
            vv += w * l * float(uc.shear_shape(xi))
        d[i], m[i], v[i] = y * 1000.0, mm, vv
    return Profile(xs, d, m, v)


def curves(inp: GeneralInput, n: int = 201) -> Profile:
    return profile(inp, np.linspace(0.0, inp.length_m, n))


def solve(inp: GeneralInput, n: int = 801) -> GeneralResult:
    l = inp.length_m
    eps = l * 1e-9
    pts, udls = _clean_loads(inp)
    extra = []
    for p in pts:
        extra += [p.pos_m, max(p.pos_m - eps, 0.0), min(p.pos_m + eps, l)]
    for u in udls:
        extra += [u.start_m, u.end_m]
    xs = np.sort(np.concatenate([np.linspace(0.0, l, n), np.array(extra, dtype=float)]))
    pr = profile(inp, xs)
    i_max = int(np.argmax(pr.deflection_mm))
    d_max = max(0.0, float(pr.deflection_mm[i_max]))
    d_at = float(xs[i_max]) if d_max > 0 else 0.0
    sag = max(0.0, float(pr.moment_nm.max()))
    hog = min(0.0, float(pr.moment_nm.min()))
    shear = float(np.abs(pr.shear_n).max())
    m_abs = max(sag, -hog)
    stress = m_abs * 1000.0 * inp.c_mm / inp.i_mm4 if inp.c_mm is not None else None
    allow = inp.strength_mpa / inp.safety_factor if inp.strength_mpa is not None else None
    util = stress / allow if stress is not None and allow is not None else None
    tau = inp.shear_factor * shear / inp.shear_area_mm2 if inp.shear_area_mm2 else None
    tau_allow = inp.shear_strength_mpa / inp.safety_factor if inp.shear_strength_mpa is not None else None
    tau_util = tau / tau_allow if tau is not None and tau_allow is not None else None
    ratio = l * 1000.0 / d_max if d_max > 0 else math.inf
    return GeneralResult(d_max, d_at, sag, hog, m_abs, shear, stress, allow, util, tau, tau_allow, tau_util,
                         ratio, ratio >= inp.limit_denominator)


# --------------------------------------------------------------------------
# Upgrades: section selection, influence line, text report
# --------------------------------------------------------------------------
def select_ipe(supports: str, length_m: float, mat: Material, point_loads=(), udl_n_per_m: float = 0.0,
               limit_denominator: int = 360, safety_factor: float = 1.5, with_self_weight: bool = True,
               n: int = 201):
    """Check every IPE size (light to heavy). Returns a list of dicts; 'ok' = deflection and bending both pass."""
    rows = []
    for name, h, iy, wy, mass, tw in IPE:
        s = ipe(name)
        g = GeneralInput(supports, length_m, mat.e_gpa, s.i_mm4, list(point_loads), udl_n_per_m,
                         self_weight_n_per_m=self_weight(mass) if with_self_weight else 0.0, c_mm=s.c_mm,
                         limit_denominator=limit_denominator, strength_mpa=mat.strength_mpa,
                         safety_factor=safety_factor, shear_area_mm2=s.shear_area_mm2,
                         shear_strength_mpa=mat.shear_strength_mpa)
        r = solve(g, n=n)
        ok = r.deflection_ok and (r.utilisation is None or r.utilisation <= 1.0) \
            and (r.shear_utilisation is None or r.shear_utilisation <= 1.0)
        rows.append({"section": name, "mass_kg_per_m": mass, "deflection_mm": r.deflection_mm,
                     "span_ratio": r.span_ratio, "utilisation": r.utilisation,
                     "shear_utilisation": r.shear_utilisation, "ok": ok})
    return rows


def lightest_ipe(*args, **kwargs):
    """Name of the lightest IPE that passes all checks, or None."""
    for row in select_ipe(*args, **kwargs):
        if row["ok"]:
            return row["section"]
    return None


def influence_deflection(supports: str, xs, positions):
    """Deflection of a unit-load beam (L = 1, EI = 1): rows = load positions, columns = xs."""
    return np.array([[_point_ei_y(supports, 1.0, 1.0, a, x) for x in xs] for a in positions])


def report_text(inp: GeneralInput, res: GeneralResult, title: str = "Beam report") -> str:
    def f(v, d=2):
        return "-" if v is None else f"{v:.{d}f}"
    return "\n".join([
        f"# {title}", "",
        f"- Supports: {inp.supports}, L = {inp.length_m} m, E = {inp.e_gpa} GPa, I = {inp.i_mm4:.4g} mm4",
        f"- Max deflection: {res.deflection_mm:.2f} mm at x = {res.deflection_at_m:.2f} m (L/{res.span_ratio:.0f}; "
        f"limit L/{inp.limit_denominator}: {'OK' if res.deflection_ok else 'FAIL'})",
        f"- Max moment: {res.max_moment_nm / 1e3:.2f} kNm (sag {res.sag_moment_nm / 1e3:.2f}, hog {res.hog_moment_nm / 1e3:.2f})",
        f"- Max shear: {res.shear_n / 1e3:.2f} kN",
        f"- Bending stress: {f(res.stress_mpa, 1)} MPa, allowable {f(res.allowable_mpa, 1)} MPa, utilisation {f(res.utilisation)}",
        f"- Shear utilisation: {f(res.shear_utilisation)}",
    ])


# --------------------------------------------------------------------------
# Upgrades v3: natural frequency and moving-load envelope
# --------------------------------------------------------------------------
_LAMBDA = {"simple": math.pi, "cantilever": 1.8751040687, "fixed": 4.7300407449, "propped": 3.9266023}


def natural_frequency(supports: str, length_m: float, e_gpa: float, i_mm4: float, mass_kg_per_m: float) -> float:
    """First bending natural frequency in Hz (distributed mass, Euler-Bernoulli):
    f = lambda^2 / (2 pi L^2) * sqrt(EI / m). Propped = fixed-pinned."""
    ei = e_gpa * 1e9 * i_mm4 * 1e-12
    lam = _LAMBDA[supports]
    return lam ** 2 / (2.0 * math.pi * length_m ** 2) * math.sqrt(ei / mass_kg_per_m)


def moving_load_envelope(supports: str, length_m: float, load_n: float, e_gpa: float, i_mm4: float,
                         n_pos: int = 61, n_x: int = 121):
    """One point load crossing the span. Returns dict with x, positions and envelopes
    (max sagging moment, max hogging moment, max |shear|, max deflection in mm)."""
    l = length_m
    ei = e_gpa * 1e9 * i_mm4 * 1e-12
    xs = np.linspace(0.0, l, n_x)
    pos = np.linspace(0.0, l, n_pos)
    mom = np.array([[_point_moment(supports, load_n, l, a, x) for x in xs] for a in pos])
    shr = np.array([[_point_shear(supports, load_n, l, a, x) for x in xs] for a in pos])
    dfl = np.array([[_point_ei_y(supports, load_n, l, a, x) / ei * 1000.0 for x in xs] for a in pos])
    return {"x": xs, "positions": pos, "moment": mom, "shear": shr, "deflection": dfl,
            "sag_env": np.maximum(mom.max(axis=0), 0.0), "hog_env": np.minimum(mom.min(axis=0), 0.0),
            "shear_env": np.abs(shr).max(axis=0), "defl_env": dfl.max(axis=0)}


# --------------------------------------------------------------------------
# Continuous beams (any number of supports, overhangs, clamped ends) - finite element solver
# --------------------------------------------------------------------------
@dataclass
class ContinuousInput:
    length_m: float
    supports_x: List[float]                 # pinned supports (vertical restraint), m from the left end
    e_gpa: float
    i_mm4: float
    point_loads: List[PointLoad] = field(default_factory=list)
    udl_n_per_m: float = 0.0                # over the whole length
    partial_udls: List[PartialUdl] = field(default_factory=list)
    self_weight_n_per_m: float = 0.0
    fixed_left: bool = False                # clamp (v = 0, rotation = 0) at x = 0
    fixed_right: bool = False               # clamp at x = L
    c_mm: Optional[float] = None
    strength_mpa: Optional[float] = None
    safety_factor: float = 1.5
    limit_denominator: int = 360


@dataclass
class ContinuousResult:
    x: np.ndarray
    deflection_mm: np.ndarray               # downward positive
    moment_nm: np.ndarray                   # sagging positive
    shear_n: np.ndarray
    reactions: dict                         # {x_m: upward reaction N}; clamp moment under key (x, "M")
    deflection_max_mm: float
    deflection_at_m: float
    sag_moment_nm: float
    hog_moment_nm: float
    max_moment_nm: float
    shear_max_n: float
    stress_mpa: Optional[float]
    allowable_mpa: Optional[float]
    utilisation: Optional[float]
    span_ratio: float
    deflection_ok: bool


def equal_spans(n_spans: int, span_m: float) -> List[float]:
    """Support positions of n equal spans (n + 1 supports)."""
    return [round(i * span_m, 9) for i in range(n_spans + 1)]


def solve_continuous(inp: ContinuousInput, n_sub: int = 6) -> ContinuousResult:
    l = inp.length_m
    ei = inp.e_gpa * 1e9 * inp.i_mm4 * 1e-12
    pts = [PointLoad(p.load_n, min(max(p.pos_m, 0.0), l)) for p in inp.point_loads if p.load_n != 0.0]
    udls = [PartialUdl(u.w_n_per_m, min(max(u.start_m, 0.0), l), min(max(u.end_m, 0.0), l)) for u in inp.partial_udls]
    udls = [u for u in udls if u.w_n_per_m != 0.0 and u.end_m > u.start_m]
    w_full = inp.udl_n_per_m + inp.self_weight_n_per_m

    key = [0.0, l] + [min(max(s, 0.0), l) for s in inp.supports_x] + [p.pos_m for p in pts]
    for u in udls:
        key += [u.start_m, u.end_m]
    nodes = set(round(k, 9) for k in key)
    max_le = l / 120.0
    base = sorted(nodes)
    for a, b in zip(base[:-1], base[1:]):
        k = int(math.ceil((b - a) / max_le))
        for i in range(1, k):
            nodes.add(round(a + (b - a) * i / k, 9))
    xs_n = np.array(sorted(nodes))
    nn = len(xs_n)
    idx = {round(x, 9): i for i, x in enumerate(xs_n)}
    ndof = 2 * nn
    k_g = np.zeros((ndof, ndof))
    f_g = np.zeros(ndof)

    def w_on(a, b):
        mid = 0.5 * (a + b)
        w = w_full
        for u in udls:
            if u.start_m - 1e-12 <= mid <= u.end_m + 1e-12:
                w += u.w_n_per_m
        return w

    elems = []
    for e in range(nn - 1):
        a, b = xs_n[e], xs_n[e + 1]
        le = b - a
        ke = ei / le ** 3 * np.array([[12, 6 * le, -12, 6 * le],
                                      [6 * le, 4 * le ** 2, -6 * le, 2 * le ** 2],
                                      [-12, -6 * le, 12, -6 * le],
                                      [6 * le, 2 * le ** 2, -6 * le, 4 * le ** 2]])
        w = w_on(a, b)
        fe = np.array([-w * le / 2, -w * le ** 2 / 12, -w * le / 2, w * le ** 2 / 12])
        d = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        k_g[np.ix_(d, d)] += ke
        f_g[d] += fe
        elems.append((a, b, le, ke, fe, d, w))
    for p in pts:
        f_g[2 * idx[round(p.pos_m, 9)]] -= p.load_n

    fixed = set()
    for s in inp.supports_x:
        fixed.add(2 * idx[round(min(max(s, 0.0), l), 9)])
    if inp.fixed_left:
        fixed |= {0, 1}
    if inp.fixed_right:
        fixed |= {2 * (nn - 1), 2 * (nn - 1) + 1}
    free = [i for i in range(ndof) if i not in fixed]
    if len(fixed) < 2 or (inp.fixed_left is False and inp.fixed_right is False and len(inp.supports_x) < 2):
        raise ValueError("Beam is not stable: add supports or clamp an end.")
    u = np.zeros(ndof)
    kff = k_g[np.ix_(free, free)]
    if np.linalg.cond(kff) > 1e14:
        raise ValueError("Beam is not stable (mechanism): check the supports.")
    u[free] = np.linalg.solve(kff, f_g[free])
    r_vec = k_g @ u - f_g

    xs, dd, mm, vv = [], [], [], []
    for (a, b, le, ke, fe, d, w) in elems:
        ue = u[d]
        fl = ke @ ue - fe                    # end forces on the element
        m_left, v_left = -fl[1], fl[0]       # sagging moment and shear at the left end
        for k in range(n_sub + 1):
            s = le * k / n_sub
            t = s / le
            n1, n2 = 1 - 3 * t ** 2 + 2 * t ** 3, s * (1 - 2 * t + t ** 2)
            n3, n4 = 3 * t ** 2 - 2 * t ** 3, s * (t ** 2 - t)
            vdef = n1 * ue[0] + n2 * ue[1] + n3 * ue[2] + n4 * ue[3]
            xs.append(a + s)
            dd.append(-vdef * 1000.0)
            vv.append(v_left - w * s)
            mm.append(m_left + v_left * s - w * s * s / 2.0)
    xs, dd, mm, vv = map(np.array, (xs, dd, mm, vv))

    reactions = {}
    for dof in sorted(fixed):
        node, kind = divmod(dof, 2)
        x_ = round(float(xs_n[node]), 6)
        if kind == 0:
            reactions[x_] = float(r_vec[dof])
        elif abs(r_vec[dof]) > 1e-9:
            reactions[(x_, "M")] = float(r_vec[dof])

    i_max = int(np.argmax(dd))
    d_max = max(0.0, float(dd[i_max]))
    sag, hog = max(0.0, float(mm.max())), min(0.0, float(mm.min()))
    m_abs = max(sag, -hog)
    stress = m_abs * 1000.0 * inp.c_mm / inp.i_mm4 if inp.c_mm is not None else None
    allow = inp.strength_mpa / inp.safety_factor if inp.strength_mpa is not None else None
    util = stress / allow if stress is not None and allow is not None else None
    # governing span ratio: largest deflection relative to the longest span
    sup = sorted(set([0.0, l] + [min(max(s, 0.0), l) for s in inp.supports_x]))
    longest = max(b - a for a, b in zip(sup[:-1], sup[1:]))
    ratio = longest * 1000.0 / d_max if d_max > 0 else math.inf
    return ContinuousResult(xs, dd, mm, vv, reactions, d_max, float(xs[i_max]) if d_max > 0 else 0.0, sag, hog,
                            m_abs, float(np.abs(vv).max()), stress, allow, util, ratio,
                            ratio >= inp.limit_denominator)


# --------------------------------------------------------------------------
# v5: load combinations, live-load patterns (continuous beams), shear deflection
# --------------------------------------------------------------------------
def load_combination(dead: float, live: float, gamma_g: float = 1.35, gamma_q: float = 1.5) -> float:
    """Factored design value gamma_g * G + gamma_q * Q (use 1.0 / 1.0 for serviceability)."""
    return gamma_g * dead + gamma_q * live


def _bays(inp: ContinuousInput):
    pts = sorted(set([0.0, inp.length_m] + [min(max(s, 0.0), inp.length_m) for s in inp.supports_x]))
    return list(zip(pts[:-1], pts[1:]))


def pattern_envelope(inp: ContinuousInput, dead_n_per_m: float, live_n_per_m: float,
                     gamma_g: float = 1.35, gamma_q: float = 1.5, n_sub: int = 4, max_bays: int = 8):
    """Envelope over every on/off pattern of the live uniform load on the bays between supports.

    Everything already in `inp` (udl, partial udls, point loads, self-weight) is permanent and is factored with
    gamma_g together with `dead_n_per_m`. Returns dict: x, m_max, m_min, v_max, v_min, d_max, d_min (deflection,
    downward positive, mm), r_max / r_min (reaction envelope per support), patterns (list of tuples of loaded bay
    indices), moments (matrix patterns x points), bays."""
    bays = _bays(inp)
    if len(bays) > max_bays:
        raise ValueError(f"At most {max_bays} bays (2^n patterns).")
    perm = gamma_g * (inp.udl_n_per_m + inp.self_weight_n_per_m + dead_n_per_m)
    pts = [PointLoad(p.load_n * gamma_g, p.pos_m) for p in inp.point_loads]
    pudl = [PartialUdl(u.w_n_per_m * gamma_g, u.start_m, u.end_m) for u in inp.partial_udls]
    patterns = [tuple(i for i in range(len(bays)) if mask >> i & 1) for mask in range(1, 2 ** len(bays))]
    patterns.insert(0, ())
    rows = []
    reacts = []
    for pat in patterns:
        live = [PartialUdl(gamma_q * live_n_per_m, bays[i][0], bays[i][1]) for i in pat] if live_n_per_m else []
        c = dataclasses.replace(inp, udl_n_per_m=perm, self_weight_n_per_m=0.0, point_loads=pts,
                                partial_udls=pudl + live)
        r = solve_continuous(c, n_sub=n_sub)
        rows.append(r)
        reacts.append({k: v for k, v in r.reactions.items() if not isinstance(k, tuple)})
    mom = np.array([r.moment_nm for r in rows])
    shr = np.array([r.shear_n for r in rows])
    dfl = np.array([r.deflection_mm for r in rows])
    keys = sorted(reacts[0].keys())
    rmat = np.array([[rc[k] for k in keys] for rc in reacts])
    return {"x": rows[0].x, "m_max": mom.max(axis=0), "m_min": mom.min(axis=0), "v_max": shr.max(axis=0),
            "v_min": shr.min(axis=0), "d_max": dfl.max(axis=0), "d_min": dfl.min(axis=0),
            "r_keys": keys, "r_max": rmat.max(axis=0), "r_min": rmat.min(axis=0),
            "patterns": patterns, "moments": mom, "bays": bays}


def shear_deflection_mm(supports: str, uniform: bool, load: float, length_m: float, e_gpa: float,
                        shear_area_mm2: float, nu: float = 0.3) -> Optional[float]:
    """Extra (Timoshenko) deflection from shear for the standard cases, mm. shear_area_mm2 = kappa * A
    (rectangle: 5/6 b h; I-beam: web area h*tw). Returns None for the propped cantilever."""
    g = e_gpa * 1e9 / (2.0 * (1.0 + nu))
    kga = g * shear_area_mm2 * 1e-6
    l = length_m
    if supports in ("simple", "fixed"):
        d = load * l * l / (8.0 * kga) if uniform else load * l / (4.0 * kga)
    elif supports == "cantilever":
        d = load * l * l / (2.0 * kga) if uniform else load * l / kga
    else:
        return None
    return d * 1000.0
