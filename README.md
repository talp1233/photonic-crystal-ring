# Squeezing in a SiN microring — Lustig lab (Technion)

Shared project repo for **Tal & Yaniv** (supervisor: Dr. Eran Lustig).

**Goal:** an on-chip source of non-degenerate two-mode squeezing in a silicon-nitride microring
(pump 1550 nm, signal 1540 nm, idler 1560 nm) — the most squeezing for the least pump power.
Eran's requirement: pump critically coupled, squeezed light overcoupled (escape efficiency ≈ 90%).
Baseline ring: Vaidya et al., *Science Advances* 2020 (Xanadu).

## Layout

| Folder | Contents |
|---|---|
| `code/` | Simulation code: `squeezing_scenarios.py` (main engine — Yaniv's `ring_three_configs.py` + new layouts + two-point/MZI coupler), MPB/MEEP scripts (`xsec.py`, `band3d.py`, `pert.py`, `defect_cavity_checks.py`) and their saved `.json` results |
| `docs/` | Write-ups: `HANDOFF_for_Yaniv.md` (current state, Nov 17), `PhC_mirror_scenarios.md` (PhC ring + mirrors, fair comparison) |
| `meetings/` | Material for / from meetings with Eran |

## Quick start

```bash
pip install numpy scipy matplotlib
python3 code/squeezing_scenarios.py --eta        # seconds
python3 code/squeezing_scenarios.py --mzi        # two-point coupler, ~2–3 min
python3 code/squeezing_scenarios.py --fair --quick
```

MPB/MEEP scripts need `pymeep` (Linux/macOS or WSL2 — see `docs/HANDOFF_for_Yaniv.md` §8):

```bash
conda create -n mp -c conda-forge python=3.11 pymeep numpy matplotlib
```

## Working together

- Pull before you start: `git pull`
- Commit small, with a clear message: `git commit -m "mzi: add Kerr shift to pump solver"`
- For bigger changes, open a branch (`git switch -c kerr`) and a pull request so the other person can look first.
- Don't commit large generated output (plots, `.h5` fields) — they're ignored by `.gitignore`. Keep small result `.json` files.
