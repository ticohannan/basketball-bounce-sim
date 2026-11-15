# User Guide — Plan

**Status:** plan (v0.2, 2026-10-20). The guide will ship with software release v1.3.0.
**Location when written:** `docs/USER_GUIDE.md`, also published on GitHub Pages.

## 1. Goals

* An undergraduate with **no prior impact-mechanics background** can:
  * install the tool,
  * run a drop,
  * interpret every plot and readout,
  * export data, and
  * make one model change,

  all within about 1 hour.
* Instructors can use the guide as a ready-made lab handout.
* The guide complements, but does not repeat, [`DESIGN.md`](DESIGN.md), which is the
  engineering reference for equations and assumptions.

## 2. Audiences

| Audience | Needs | Primary sections |
|---|---|---|
| Student user | Run it, read the plots, export data | 3–6, 8 |
| Student modifier | Change ranges or physics safely | 7, 9 |
| Instructor | Learning objectives, exercises, expected answers | 8, Appendix B |
| Contributor | Test and release process | Link to DESIGN §12 |

## 3. Proposed outline

1. **Introduction.** What the tool does, and what it doesn't (from DESIGN §1 and §6), in
   plain language.
2. **Installation.**
   * 2.1 WinPython 3.12 on Windows (step-by-step with screenshots).
   * 2.2 Anaconda / pip on Windows, macOS and Linux.
   * 2.3 JupyterLab vs classic Notebook vs VS Code. Widget support notes: if the sliders
     don't appear, check that ipywidgets is installed and trust the notebook.
3. **First run (tutorial).**
   * Reproduce the 12 ft / 8 psi example.
   * Annotated screenshot of the panel, with numbered call-outs for every control.
4. **Controls reference.** One subsection per slider or toggle:
   * the physical meaning,
   * the published range and its source,
   * what to expect when it is changed.

   This includes:
   * the initial-condition toggle ("Drop from rest" / "Initial velocity"), the sign
     convention (+ up) and the m/s / ft/s toggle;
   * the ft/m and psi/kPa toggles;
   * the underside-of-ball height convention;
   * the air-drag switch.
5. **Reading the outputs.**
   * 5.1 Readout panel, line by line.
   * 5.2 Plots (a)–(f): what each axis means and what features to look for. Examples:
     the peak force comes *before* maximum deflection; the area of the hysteresis loop is
     the energy lost.
   * 5.3 Rule-check verdicts, and why FIBA 2026 measures the rebound to the underside of
     the ball while FIBA 2018 and NFHS measure to the top.
6. **Exporting and using data.**
   * CSV structure and the metadata lines.
   * Opening the file in Excel, and in pandas with `comment="#"`.
   * Plotting the exported data; comparing it with force-plate or piezo measurements.
7. **Scripting without the GUI.**
   * `simulate_drop`, `sweep_ranges`, `calibrate_damping`.
   * A worked parameter study: contact time vs pressure.
8. **Guided exercises** (with expected results in Appendix B):
   1. Energy budget at maximum deflection. Why is the shell share so small?
   2. Under- vs over-inflation (5 vs 12 psi): contact time, peak force, COR.
   3. Does rubber hardness matter? (The storage-number argument.)
   4. Drop-height sweep: why contact time barely changes while peak force grows roughly
      with √h.
   5. Turn damping off: check that COR = 1 and the energy balance closes.
   6. Re-calibrate to the NFHS test instead of FIBA. Compare the λ values.
   7. Lab tie-in: measure a real ball's rebound (as in [13]) and fit λ.
   8. Dribble vs drop: throw the ball down at 3–6 m/s from waist height. Compare the
      equivalent drop height and the peak force with a 12 ft drop. Then explain why an
      upward throw with drag hits more softly than a downward throw at the same speed.
9. **Modifying the model.** Recipes, each with a "verify with pytest" step:
   * different ball size,
   * compliant floor,
   * Kelvin–Voigt damping,
   * isothermal gas.
10. **Troubleshooting and FAQ.**
    * Widgets show as text.
    * Slow first run (the sweep).
    * CSV saved "somewhere else" (the working directory).
    * The "Ball never left the floor" error (damping too high).
11. **Glossary.** Gauge vs absolute pressure, COR, Shore A, adiabatic, hysteresis, impulse.

* **Appendix A** — Units and conversions.
* **Appendix B** — Instructor answer key.
* **Appendix C** — References (link to [`REFERENCES.md`](REFERENCES.md)).

## 4. Assets to produce

| Asset | Source | Status |
|---|---|---|
| Panel screenshots: wide, narrow, metric, initial velocity | `docs/images/` | ✅ done (v1.2.0) |
| Example bounce and sweep figures | `docs/images/example_*.png` | ✅ done |
| Example CSV | `examples/` | ✅ done |
| Annotated panel screenshot with numbered call-outs | New | ☐ |
| Installation screenshots (WinPython) | New, captured on Windows | ☐ |
| Exercise notebooks (`examples/exercise_*.ipynb`) | New | ☐ |
| Short GIF of slider interaction | New | ☐ |

## 5. Style and conventions

* Second person, short steps, one action per numbered step.
* SI units first, with imperial units in parentheses (matching the panel).
* Every number quoted in the guide must be reproducible with the stated settings and
  software version.
* Use the same reference numbering as `REFERENCES.md` (ASME style).

## 6. Review and maintenance

* Technical review: Tico.
* Usability review: one undergraduate who has never seen the tool completes Sections 2–6
  unaided. Record where they get stuck.
* Re-check all screenshots and quoted numbers whenever the MINOR or MAJOR version changes.
  The release checklist is in DESIGN §12.

## 7. Milestones

| Milestone | Content | Target version |
|---|---|---|
| M1 | Sections 1–6 + troubleshooting | v1.3.0 |
| M2 | Exercises 8.1–8.6 and 8.8 + answer key | v1.4.0 |
| M3 | Modification recipes + lab tie-in exercise | v1.5.0 |
