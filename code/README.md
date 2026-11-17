# code/

Simulation files (run from the repo root, e.g. `python3 code/squeezing_scenarios.py --eta`):

- `squeezing_scenarios.py` — main engine (map at top: PART 1 … PART 11; all layouts in `LAYOUT_SPEC`, PART 7). Needs numpy + scipy.
- `defect_cavity_checks.py` — subcommands `gap` / `cmt` / `cavity` / `field` / `overlap` / `summary`. Reads and writes its JSON results next to itself (`code/`).
  - `check1.json`, `cav_abrupt.json`, `cav_taper.json`, `overlap.json` — saved results read by `summary` (`cmt.json` is not included; regenerate with `cmt`).
- `xsec.py` — MPB waveguide cross-section (needs meep/mpb)
- `band3d.py` — MPB 3D corrugated unit cell (needs meep/mpb)
- `pert.py` — MPB + vector perturbation theory (needs meep/mpb)
- `ring_three_configs.py` — Yaniv's original engine (for reference)
