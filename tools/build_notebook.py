"""
Regenerate basketball_bounce_sim.ipynb from basketball_bounce_sim.py.

The .py file is the single source of truth.  The notebook is a convenience
copy: one Markdown cell of instructions + one code cell containing the whole
module, so it runs on its own (no %run or import needed).

Usage (from the repository root):
    python tools/build_notebook.py
"""
import os
import re

import nbformat as nbf

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "basketball_bounce_sim.py")
DST = os.path.join(ROOT, "basketball_bounce_sim.ipynb")

code = open(SRC, encoding="utf-8").read()
version = re.search(r'__version__ = "([^"]+)"', code).group(1)

intro = f"""# Basketball drop-impact simulator (v{version})

Run the code cell below (**Shift+Enter**). An interactive panel appears with sliders for
drop height (ft or m), inflation pressure (psi or kPa), ball mass and size, rubber
hardness, shell thickness and damping.

* Plots redraw when you **release** a slider. Current values are printed under each slider.
* **Calibrate λ** sets the damping so the ball passes the FIBA 2026 rebound test
  (1800 mm drop → 1035–1085 mm rebound, measured to the underside of the ball).
* **Export CSV** writes the time history (and optionally the min/max range sweep) next to this notebook.

The first run takes a few seconds while the regulation-range sweep is computed.
Equations, assumptions and references: see `docs/DESIGN.md` and `docs/REFERENCES.md`.

*This notebook is generated from `basketball_bounce_sim.py` by `tools/build_notebook.py`;
edit the .py file, then rebuild.*"""

nb = nbf.v4.new_notebook()
nb.cells = [nbf.v4.new_markdown_cell(intro), nbf.v4.new_code_cell(code)]
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nb.metadata["language_info"] = {"name": "python"}
nbf.write(nb, DST)
print(f"wrote {DST} (v{version})")
