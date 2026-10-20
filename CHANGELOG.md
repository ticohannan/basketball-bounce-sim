# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] — 2026-10-20

### Added
- **Initial-velocity option.** A toggle switches between "Drop from rest" and
  "Initial velocity". The new starting-velocity slider runs from −15 to +15 m/s (+ = up),
  shown in **m/s or ft/s**.
- Exact quadratic-drag solutions for thrown-up and thrown-down starts:
  `impact_speed(h, ball, v_init)` and the new helper `rise_height()`.
- New readouts: start condition, apex height and equivalent drop height. The same values
  are written to the CSV metadata.
- The range sweep and its plots use the chosen starting velocity.
- 10 new tests (23 in total).
- `LICENSE` (BSD-3-Clause). Author and repository details added to `CITATION.cff`, the
  README and the code header.
- Screenshot of initial-velocity mode.

### Changed
- `params` accepts the optional key `v_init_m_s`. It defaults to 0, so existing scripts
  behave exactly as before.
- Removed the verification log from `docs/REFERENCES.md`.

## [1.1.0] — 2026-10-20

### Added
- Unit toggles: drop height in **ft or m**, inflation pressure in **psi or kPa**. Values are
  converted, so the physics is unaffected by the display unit.
- Slider values are printed **below** each slider in both unit systems, so they stay readable
  in narrow windows. The control columns now wrap on small screens.
- `__version__` constant. The version is shown in the panel title and written to the CSV
  metadata.
- FIBA 2026 equipment rule (valid from 1 Oct 2026): rebound 1035–1085 mm measured to the
  **underside** of the ball. This is now the default rule for damping calibration.
- FIBA 2018 test kept as a legacy rule check.
- Pressure in kPa added to the CSV metadata. The range-sweep CSV now includes heights in
  metres.
- Verification test suite (`tests/test_physics.py`, 13 tests).
- Documentation: design document, ASME-format references with verified DOIs, user guide
  plan, GitHub Pages site with screenshots.
- `CITATION.cff`, `requirements.txt`, and `tools/build_notebook.py` (generates the notebook
  from the `.py` file).

### Changed
- The NFHS rule test is now labelled "NFHS". The earlier claim that NCAA uses the identical
  test was removed because it could not be confirmed.
- Gent (1958) is now cited by its Crossref-registered *Rubber Chem. Technol.* record.

## [1.0.0] — 2026-10-20

### Added
- First release: adiabatic gas-spring and thin-shell contact model with Hunt–Crossley damping.
- Air drag during fall and rebound.
- Damping calibrated to the FIBA 2018 rebound test.
- Interactive ipywidgets panel, four-panel bounce figure and regulation-range envelope.
- CSV export.
