# BeamDeflection - Jupyter Notebook

Beam deflection, bending moment, shear force and stress checks for single-span beams. This is the Python / Jupyter version of the **BeamDeflection Android app** (Kotlin + Jetpack Compose), as a **flat project: no sub-folders**.

![Deformed beam in 3D](fig_3d_deformed_beam.png)

## Features

- Continuous beams (any spans, overhangs, clamped ends) and a Streamlit dashboard
- Automatic lightest-IPE selection, influence surfaces, text report
- Vibration frequency, moving loads, sensitivity, Monte Carlo, Pareto front

- 4 support types: simply supported, cantilever, fixed-fixed, propped cantilever
- Point loads anywhere on the span, full or partial uniform loads, self-weight (superposition)
- 14 materials and the IPE 100-300 catalogue; rectangle, round bar, box and tube sections
- Deflection, span/deflection limit (L/180 to L/500), bending moment, shear force, bending and shear utilisation
- 2D diagrams, **3D diagrams, 3D charts, 3D graphs** and **comparison tables** in the notebook

## Files (all in the project root)

| File | Purpose |
|---|---|
| `BeamDeflection.ipynb` | Main notebook: example, diagrams, 3D charts, tables, checks |
| `beam_core.py` | Solver, materials, sections (port of the Kotlin domain and data code), continuous-beam solver |
| `beam_plots.py` | Shared matplotlib figures (profiles, 3D beam, reactions) |
| `app.py` | Streamlit dashboard (`streamlit run app.py`) |
| `test_beam_core.py` | Tests against closed-form formulas |
| `fig_*.png` | Figures written by the notebook |
| `example_report.md`, `example_report.pdf` | Sample text and PDF reports |
| `example_profile.csv` | Deflection / moment / shear profile of the example beam |
| `requirements.txt` | numpy, matplotlib, pandas, jupyter, streamlit |

## Quick start

```bash
pip install -r requirements.txt
jupyter notebook BeamDeflection.ipynb
streamlit run app.py
python test_beam_core.py
```

```python
import beam_core as bc
sec, mat = bc.ipe("IPE 200"), bc.material("S275")
beam = bc.GeneralInput("simple", 5.0, mat.e_gpa, sec.i_mm4, c_mm=sec.c_mm,
                       point_loads=[bc.PointLoad(20_000, 2.0)], udl_n_per_m=6_000,
                       self_weight_n_per_m=bc.self_weight(sec.mass_kg_per_m),
                       strength_mpa=mat.strength_mpa)
print(bc.solve(beam))
```

## Diagrams

### Workflow

```mermaid
flowchart LR
    A[Inputs: supports, L, E, section, loads] --> B[beam_core.GeneralInput]
    M[Material library] --> A
    S[Section library / IPE] --> A
    B --> C[Superposition solver]
    C --> D[Deflection, moment, shear profiles]
    D --> E[Checks: L/n limit, bending, shear]
    D --> F[2D and 3D plots]
    E --> G[Result tables]
```

### Module structure

```mermaid
classDiagram
    class GeneralInput {
        supports
        length_m
        e_gpa
        i_mm4
        point_loads
        udl_n_per_m
    }
    class LoadCase {
        shape()
        moment_shape()
        shear_shape()
    }
    class Material {
        e_gpa
        density
        strength_mpa
    }
    class SectionProps {
        i_mm4
        c_mm
        shear_area_mm2
    }
    class GeneralResult {
        deflection_mm
        max_moment_nm
        utilisation
    }
    GeneralInput --> LoadCase : self-weight and UDL
    GeneralInput ..> Material
    GeneralInput ..> SectionProps
    GeneralInput --> GeneralResult : solve()
```

### Support types

```
Simply supported      Cantilever           Fixed-fixed          Propped cantilever
 v P                   v P                  v P                  v P
 /\======o            ##========           ##========##         ##========o
```

### 2D diagrams (example beam)

![Deflection, moment and shear](fig_2d_profiles.png)

## 3D charts and graphs

**Deflection vs span and load**

![3D surface](fig_3d_surface_span_load.png)

**Support types compared (same total load)**

![3D support comparison](fig_3d_support_comparison.png)

**IPE size x material**

![3D IPE x material](fig_3d_ipe_material.png)

**Bending utilisation vs span and depth**

![3D utilisation](fig_3d_utilisation.png)

**Moment-shear-position curve**

![3D curve](fig_3d_moment_shear_curve.png)

**Material space (density x stiffness x strength)**

![3D materials](fig_3d_material_ashby.png)

**Influence surfaces, all four supports**

![3D influence](fig_3d_influence_surfaces.png)

**Bending moment ribbons**

![3D ribbons](fig_3d_moment_ribbons.png)

**Lightest passing IPE by span and load**

![3D IPE selection](fig_3d_ipe_selection.png)

## Upgrades in this version

- `select_ipe()` / `lightest_ipe()`: checks every IPE size for deflection limit, bending and shear (self-weight included)
- `influence_deflection()`: deflection influence surfaces for any support type
- `report_text()`: short markdown report (the notebook writes `example_report.md`)
- 4 new 3D charts, more tests
- v5: `pattern_envelope()` (live-load patterns, ULS/SLS), `shear_deflection_mm()`, `make_pdf_report()`, dashboard tab and PDF download, 3 more tests
- v4: `solve_continuous()`, `app.py` dashboard, `beam_plots.py`, 3 more tests
- v3: `natural_frequency()` (first bending mode, 4 supports), `moving_load_envelope()` (moment, shear, deflection envelopes)
- v3 notebook sections: sensitivity (tornado), moving load, vibration, Monte Carlo, mass-vs-deflection Pareto front, CSV export (`example_profile.csv`)

**Moving-load moment surface**

![3D moving load](fig_3d_moving_load_moment.png)

**Natural frequency vs span and IPE size**

![3D frequency](fig_3d_natural_frequency.png)

## Continuous beams (v4)

`solve_continuous()` is a finite-element solver for beams on any number of supports: several spans, overhangs, clamped ends, point loads, full and partial uniform loads, self-weight, support reactions. Checked against closed-form results (two-span beam: middle reaction 1.25 wL, hogging moment wL²/8).

![Continuous profiles](fig_continuous_profiles.png)

![3D continuous beam](fig_3d_continuous_beam.png)

**One span vs several spans (same total length)**

![3D spans](fig_3d_continuous_comparison.png)

**Best interior support positions (3-span beam)**

![3D support optimisation](fig_3d_support_optimisation.png)

## Load combinations and live-load patterns (v5)

`pattern_envelope()` places the live load on every combination of spans and returns min/max envelopes of moment, shear, deflection and support reactions (ULS 1.35 G + 1.5 Q, SLS 1.0 / 1.0; factors are inputs). Checked on a two-span beam: hogging wL²/8 with both spans loaded, sagging 49/512 wL² with one span loaded.

![Envelopes](fig_pattern_envelope.png)

![3D patterns](fig_3d_pattern_moment.png)

## Shear deformation (v5)

`shear_deflection_mm()` adds the Timoshenko shear term for the standard cases. It matters for short, deep beams.

![3D shear share](fig_3d_shear_share.png)

## PDF report (v5)

`beam_plots.make_pdf_report()` writes an A4 PDF (inputs, results, diagrams, 3D beam). The notebook writes `example_report.pdf`; the dashboard has a PDF download button.

## Streamlit dashboard (v4)

```bash
pip install -r requirements.txt
streamlit run app.py
```

```mermaid
flowchart TD
    U[Sidebar: material, section, safety, limit] --> A[app.py]
    L[Main: spans, supports, loads tables] --> A
    A -->|single span| S[bc.solve / bc.curves]
    A -->|continuous| C[bc.solve_continuous]
    S --> R[Results, diagrams, 3D view]
    C --> R
    R --> T[Tools: material comparison, lightest IPE, frequency, CSV and report download]
```

Features: single-span and continuous beams, a **Load patterns** tab (ULS / SLS envelopes), PDF report, IPE / rectangle / round / tube / custom sections, editable load tables, utilisation checks, diagrams, 3D view, material comparison, lightest-IPE finder, natural frequency, CSV and markdown download.

## Analysis charts (v3)

**Sensitivity (tornado)**

![Tornado](fig_tornado_sensitivity.png)

**Moving-load envelopes**

![Envelopes](fig_moving_load_envelope.png)

**Monte Carlo uncertainty**

![Monte Carlo](fig_monte_carlo.png)

**Mass vs deflection, Pareto front**

![Pareto](fig_pareto_mass_deflection.png)

## Comparison tables

### Example result (IPE 200, S275, L = 5 m)

| Quantity                 |    Value |
|--------------------------|----------|
| Max deflection (mm)      |    24.48 |
| at x (m)                 |     2.43 |
| Span / deflection        |   204    |
| Limit L/360 OK           | False    |
| Max sagging moment (kNm) |    42.66 |
| Max shear (kN)           |    27.55 |
| Bending stress (MPa)     |   219.6  |
| Allowable (MPa)          |   183.3  |
| Bending utilisation      |     1.2  |
| Shear utilisation        |     0.22 |

### Standard load cases

| Case                                      | Max deflection        |   d coeff |   M coeff |   V coeff | d vs simple UDL   |
|-------------------------------------------|-----------------------|-----------|-----------|-----------|-------------------|
| Simply supported - point load at centre   | d = PL^3 / 48EI       |   0.02083 |    0.25   |    0.5    | -                 |
| Simply supported - uniform load           | d = 5wL^4 / 384EI     |   0.01302 |    0.125  |    0.5    | 1.0               |
| Cantilever - point load at free end       | d = PL^3 / 3EI        |   0.33333 |    1      |    1      | -                 |
| Cantilever - uniform load                 | d = wL^4 / 8EI        |   0.125   |    0.5    |    1      | 9.6               |
| Fixed both ends - point load at centre    | d = PL^3 / 192EI      |   0.00521 |    0.125  |    0.5    | -                 |
| Fixed both ends - uniform load            | d = wL^4 / 384EI      |   0.0026  |    0.0833 |    0.5    | 0.2               |
| Propped cantilever - point load at centre | d ~ 0.00932 PL^3 / EI |   0.00932 |    0.1875 |    0.6875 | -                 |
| Propped cantilever - uniform load         | d ~ 0.00542 wL^4 / EI |   0.00542 |    0.125  |    0.625  | 0.42              |

### Support types, 20 kN total load (IPE 200, S275, L = 4 m)

| Support    |   δ point (mm) |   δ UDL (mm) |   M point (N·m) |   M UDL (N·m) |
|------------|----------------|--------------|-----------------|---------------|
| simple     |           6.54 |         4.08 |           20000 |      10000    |
| cantilever |         104.57 |        39.21 |           80000 |      40000    |
| fixed      |           1.63 |         0.82 |           10000 |       6666.67 |
| propped    |           2.92 |         1.7  |           15000 |      10000    |

### Materials

| Material                        |   E (GPa) |   Density (kg/m3) | Strength (MPa)   | Shear strength (MPa)   |   Specific stiffness E/rho (MNm/kg) |
|---------------------------------|-----------|-------------------|------------------|------------------------|-------------------------------------|
| Structural steel S235           |     210   |              7850 | 235.0            | 141.0                  |                                26.8 |
| Structural steel S275           |     210   |              7850 | 275.0            | 165.0                  |                                26.8 |
| Structural steel S355           |     210   |              7850 | 355.0            | 213.0                  |                                26.8 |
| Stainless steel 304             |     193   |              8000 | 215.0            | 129.0                  |                                24.1 |
| Aluminium 6061-T6               |      68.9 |              2700 | 276.0            | 166.0                  |                                25.5 |
| Aluminium 6063-T5               |      68.9 |              2700 | 145.0            | 87.0                   |                                25.5 |
| Titanium Ti-6Al-4V              |     113.8 |              4430 | 880.0            | 528.0                  |                                25.7 |
| Copper (annealed)               |     117   |              8960 | 70.0             | 42.0                   |                                13.1 |
| Grey cast iron                  |     110   |              7200 | -                | -                      |                                15.3 |
| Concrete C30/37                 |      33   |              2400 | -                | -                      |                                13.8 |
| Timber C24 (softwood)           |      11   |               420 | 24.0             | 4.0                    |                                26.2 |
| Timber D40 (hardwood)           |      11   |               660 | 40.0             | 4.0                    |                                16.7 |
| Glulam GL24h                    |      11.5 |               420 | 24.0             | 3.5                    |                                27.4 |
| CFRP laminate (quasi-isotropic) |      70   |              1600 | -                | -                      |                                43.8 |

### IPE sections (S275, simply supported, L = 4 m, 5 kN/m + self-weight)

| Section   |   h (mm) |   I (cm4) |   Mass (kg/m) |   Deflection (mm) |   L/deflection |   Stress (MPa) |   Utilisation |
|-----------|----------|-----------|---------------|-------------------|----------------|----------------|---------------|
| IPE 100   |      100 |       171 |           8.1 |             47.15 |             85 |          297   |          1.62 |
| IPE 120   |      120 |       318 |          10.4 |             25.47 |            157 |          192.5 |          1.05 |
| IPE 140   |      140 |       541 |          12.9 |             15.04 |            266 |          132.7 |          0.72 |
| IPE 160   |      160 |       869 |          15.8 |              9.42 |            425 |           94.9 |          0.52 |
| IPE 180   |      180 |      1317 |          18.8 |              6.25 |            640 |           70.9 |          0.39 |
| IPE 200   |      200 |      1943 |          22.4 |              4.26 |            938 |           53.7 |          0.29 |
| IPE 220   |      220 |      2772 |          26.2 |              3.01 |           1329 |           41.7 |          0.23 |
| IPE 240   |      240 |      3892 |          30.7 |              2.16 |           1850 |           32.7 |          0.18 |
| IPE 270   |      270 |      5790 |          36.1 |              1.47 |           2725 |           25   |          0.14 |
| IPE 300   |      300 |      8356 |          42.2 |              1.03 |           3889 |           19.4 |          0.11 |

Table images: `fig_table_load_cases.png`, `fig_table_materials.png`, `fig_table_ipe.png`.

## Assumptions and limits

Linear-elastic Euler-Bernoulli beam, constant section, small deflections, single span, strong-axis bending. Material and section values are typical handbook values. Shear is a simplified check; buckling and connections are not checked. Not a substitute for a code-compliant structural design.

## License

MIT - see [LICENSE](LICENSE). © 2026 Harpreet Singh
