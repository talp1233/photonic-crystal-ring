# Handoff: squeezing in a SiN microring — where Tal's side of the project stands (Nov 15, 2026; §6.2 added Nov 17)

*Written for Yaniv and for the Claude assistant helping him. It explains the files Tal is sending, what was found, and the direction being checked now. Everything numerical below comes from Yaniv's own engine (`ring_three_configs.py`) or from MPB runs; nothing was changed in the engine's physics. The two-point coupler model (PART 9b) is new code alongside it.*

---

## 1. One-paragraph summary

The five bus-reflector layouts on a **regular** ring (Yaniv's Nov 8 study) were extended with three **PhC ring + reflector** layouts, as Eran asked. Result: **a PhC ring with a single pump reflector on the right is the first layout that beats the plain ring in the fair comparison** (intrinsic Q fixed, each layout at its best η), but only at high targets: 0.61× the plain ring's pump at 10 dB and 0.40× at 12 dB with ideal reflectors (0.87× / 0.75× with 2% reflector loss). Below about 8 dB it loses (1.3–1.8×). Adding a full reflector or a second (squeeze) reflector does not help. An MPB check of the real Vaidya waveguide confirms the engine's inputs (n_g, γ) and shows that a uniform PhC corrugation does **not** shrink the mode volume, so Eran's "better confinement" idea needs a different PhC design (a localized/defect mode), which is the open direction now. **Defect-cavity checks (§6.1) say that route does not pay as designed today:** the mode volume does shrink about 3× (F ≈ 3), but radiation loss at the defect and the frequency mismatch eat the gain; it would need Q_rad ≳ 4×10⁶ on all three modes and pump/signal/idler matched within ~1 GHz. **The two-point (MZI-coupled) coupler is now a full model in the engine (§6.2):** it beats the plain ring at every target and stays robust with realistic imperfections (1% bypass loss, coupler dispersion, 10° delay and heater errors): **0.89× at 4 dB, 0.53× at 8 dB, 0.38× at 10 dB, 0.26× at 12 dB**, vs 1.62× / 0.91× / 0.61× / 0.40× for the PhC ring + pump reflector with ideal reflectors.

---

## 2. Files Tal is sending

| File | What it is | Run |
|---|---|---|
| `squeezing_scenarios.py` | **Yaniv's engine merged with the new layouts and the comparison tools.** Single file, no other file needed. Pump solver, Bogoliubov squeezer and `run()` are Yaniv's code, unchanged in physics. | `python3 squeezing_scenarios.py --eta` (seconds) · `--power` (~1–2 min) · `--fair` (~30 min) · `--fair --quick` · `--mzi` (two-point coupler, §6.2, ~2–3 min) · `--loss 0.005` |
| `xsec.py` | MPB: cross-section of the Vaidya waveguide (1000 × 800 nm SiN in SiO₂) → n_eff, n_g, A_eff, γ, and grating splitting per nm (scalar estimate). | needs MPB (see §8) |
| `band3d.py` | MPB: direct 3D unit cell of a sinusoidally corrugated waveguide; gap at the Bragg point. | `python band3d.py <amplitude_um> <res_yz> <res_x> <slices>`, e.g. `0.01 48 240 64` |
| `pert.py` | MPB + vector perturbation theory for the same splitting (includes the sign flip of the longitudinal field Ez in the standing waves). | |
| `defect_cavity_checks.py` | **The defect-cavity (mode-volume) checks of §6.1, one file with subcommands:** `gap` (MPB 3D band gap + convergence + period retune), `cmt` (coupled-mode defect estimate), `cavity` (MEEP 2D defect cavity → frequencies, Q_rad), `field` + `overlap` (mode profiles from MEEP → F), `summary` (all tables + required pump vs the plain ring; numpy only). Saved results: `check1.json`, `cav_abrupt.json`, `cav_taper.json`, `overlap.json`. | `python defect_cavity_checks.py summary` (seconds, reads the JSONs); other subcommands need MEEP, ~10 min–1 h each |

The file has a map at the top ("PART 1 … PART 11"). **All layouts are defined in one place: `LAYOUT_SPEC` in PART 7.**

### Layout names (old → new)

| New key | Layout | Useful port | Yaniv's old name |
|---|---|---|---|
| `plain` | Plain ring | right | `alone` (scenario 1) |
| `ring_P` | Ring + pump reflector (R) | right | `pumpR` (2) |
| `phc` | PhC ring, no reflector | right | `sim1` (3) |
| `ring_P_S` | Ring + pump refl. (R) + squeeze refl. (L) | right | `sim2` (4) |
| `ring_F` | Ring + full reflector (R) | left | `sim3` (5) |
| `phc_P` | **PhC ring + pump reflector (R)** | both (take the better) | new (6) |
| `phc_F` | PhC ring + full reflector (R) | left | `aside` (7) |
| `phc_P_S` | PhC ring + pump refl. (R) + squeeze refl. (L) | right | new (8) |
| — | Two-point / MZI-coupled ring (Tal's proposal) | right | own model, `mzi_*()` in PART 9b (§6.2); the analytic bound `two_point_*()` in PART 9 is its ideal limit |

---

## 3. What was changed relative to `ring_three_configs.py`

1. **Fair comparison at fixed intrinsic Q** (8×10⁵, Vaidya), not fixed loaded Q: `set_ring(eta)` sets Q_loaded = Q_int(1−η).
2. **Bug fix needed for any η sweep:** the engine computed `half_linewidth` once at import. `set_ring()` now updates it with Q. Without this, x scales as ηQ instead of ηQ², and the η sweep is wrong.
3. **PhC splitting fixed in absolute terms** (10 Vaidya linewidths ≈ 10 GHz). In loaded linewidths it shrinks as η grows (`PHC_BETA` scales with Q_loaded).
4. **Reflector convention:** reflect band R = 1 − 1.5L, T = 0.5L, A = L (L = 2% → 97/1/2%, Yaniv's placeholder); pass band T = 1 − L. This convention reproduces Yaniv's mirror-loss table.
5. **`best_pump()`**: coarse grid over laser detuning D and right-reflector phase φ, then Nelder–Mead. Objective = maximum pump absorbed.
6. **Analysis tools added:** squeezing at the useful port with best/worst/"in-step" phase of returning pairs; threshold; least pump for a target (bisection below threshold); best η for a target (coarse 0.01 then fine 0.001); best η for a power.
7. **Two-point coupler bound:** analytic, x_bound = x_plain/(4η(1−η)): all pump absorbed, one direction, no bypass loss. **Now modelled in the engine (PART 9b, §6.2).**

**Validation done:** the merged file reproduces Yaniv's 100 mW table (all five layouts), his mirror-loss table (within ~0.03×) and all numbers below.

---

## 4. Main results

### 4.1 At Vaidya's ring (η = 0.75, 100 mW, 2% reflectors), squeezing in dB

| Layout | Left | Right | Pump absorbed | x CCW / CW |
|---|---|---|---|---|
| Plain ring | 0 | **4.41** | 75.0 mW | 0.44 / 0 |
| PhC ring + pump reflector | **3.28** | **2.89** | 95.4 mW | 0.26 / 0.30 |
| PhC ring + full reflector | **3.39** | 0.02 | 95.4 mW | 0.26 / 0.30 |
| PhC ring + pump & squeeze refl. | 0.02 | **3.00** | 93.5 mW | 0.26 / 0.29 |
| Two-point coupler (ideal bound) | 0 | **5.23** | 100 mW | 0.59 / 0 |

### 4.2 Fair comparison: pump needed relative to the plain ring at its best η

Plain ring reference: 46 / 75 / 186 / 454 / 1113 / 2748 mW for 3 / 4 / 6 / 8 / 10 / 12 dB (best η 0.62 → 0.94).

| Layout (ideal reflectors) | 4 dB | 8 dB | 10 dB | 12 dB |
|---|---|---|---|---|
| Ring + pump + squeeze refl. | 0.99× | 0.88× | 0.88× | 0.92× |
| **PhC ring + pump reflector** | 1.62× | 0.91× | **0.61×** | **0.40×** |
| PhC ring + full reflector | 1.61× | 0.91× | 0.61× | 0.40× |
| PhC + pump & squeeze refl. | 1.62× | 0.92× | 0.62× | 0.41× |
| Two-point coupler (bound) | 0.85× | 0.50× | 0.34× | 0.23× |

With 2% reflectors the PhC + pump reflector gives 0.87× at 10 dB and 0.75× at 12 dB; with 0.5% it gives 0.64× / 0.43×.

### 4.3 Why: three mechanisms

1. **The reflector phase is a separate coupling knob for the pump.** The PhC makes the pump a standing wave that leaks both ways (κ_ex/2 each). The right reflector sends the right-going leak back into the same standing wave, giving one port with κ_eff = κ_ex(1 + cos φ). Choosing φ gives critical coupling for the pump (cos φ = (1−η)/η − 1; at η = 0.75 this predicts 0.73π, the optimiser finds 0.71π). Signal and idler are not split, so their η is still set by the gap. This is "pump critical, squeezing overcoupled", which is Eran's requirement.
2. **The standing wave halves the drive per direction.** x_per_direction / x_plain = 1/(8η(1−η)), which crosses 1 at η ≈ 0.85. That is why it loses at low targets and stays ~1.7–1.8× above the two-point bound.
3. **Squeezed light sent back through the ring is erased.** A squeezer maps the quiet quadrature as V_out = t²V_in + ℓ², with t = 2η/(1+x) − 1 and ℓ² = 4η(1−η)/(1+x)². At the optimal operating points t² ≤ 0.005 while ℓ² ≈ 0.04–0.16, so the returning half contributes nothing. This is why `phc_F` and `phc_P_S` equal `phc_P`.

### 4.4 Checks that changed earlier expectations

- **Raising intrinsic Q does not rescue the cascade.** Scan at Q_int = 0.8 / 2 / 3 ×10⁶: all powers scale as κ₀² (10 dB plain ring: 1113 → 178 → 79 mW), but the best η per target does not move. So t² stays ≈ 0 and the cascade layouts never beat `phc_P`. The cascade only starts to pay at forced η ≥ 0.98, where the ring needs more pump than it saves.
- **"Pump-mirror phase in step":** for `phc_F` the returning pair phase φ_s + φ_i − 2φ_p is locked to the pump (same path), so the reflector distance cannot tune it; only the Bragg reflector's spectral phase sets it. Best vs in-step differ by ≤ 0.01× in the fair comparison.

---

## 5. MPB results (Vaidya waveguide 1000 × 800 nm, SiN n = 1.996, SiO₂ n = 1.444, 1550 nm)

| Quantity | MPB | Engine assumption |
|---|---|---|
| n_g | 2.10 | 2.09 (from FSR = 190 GHz) |
| A_eff | 0.94 µm² | — |
| γ (n₂ = 2.4×10⁻¹⁹ m²/W assumed) | 1.04 W⁻¹m⁻¹ | 1.0 (`gamma_NL = 1`) |
| Bragg period | 437 nm (~1726 periods in the R = 120 µm ring) | — |
| Splitting per nm of sinusoidal width modulation (both walls) | ~6 GHz/nm (direct 3D); ≤ ~10 GHz/nm (vector perturbation) | — |

**Meaning:** the assumed 10 GHz splitting corresponds to ~1–1.7 nm corrugation per wall (~2–3.5 nm if only the inner wall is corrugated). A uniform 1–2 nm corrugation leaves A_eff and the mode length unchanged, so **the mode volume does not shrink (F ≈ 1)** in the PhC we modelled. The 6 vs 10 GHz/nm gap between methods is not resolved (3D resolution vs wall sampling). No ring curvature was included.

---

## 6. The direction being checked now

Eran said a PhC ring localizes the mode with better confinement, and asked to think about the mode volume, the PhC size and MPB. In the model, x ∝ U/(V_eff·κ), so a mode volume smaller by a factor F divides every required pump by F. For example, F = 1.6 would bring the PhC + pump reflector to parity with the plain ring already at 4 dB.

But:
1. **Only localizing the pump does not help.** The coupling scales with the overlap ∫u_p²u_s*u_i* dV. If the pump is localized but signal and idler are spread over the ring, the overlap stays ~1/L. All three modes must share the small volume.
2. **A uniform corrugation does not localize** (MPB: F ≈ 1).
3. A real gain needs a **defect (localized) PhC mode** or a strongly modulated region, and that changes Q, the bus coupling and the frequency matching.

### 6.1 Defect-cavity checks (Nov 15) — results

Design checked: the Vaidya waveguide with a **strong** sidewall grating (±200 nm per wall, period ≈ 440 nm) on both sides of a plain **defect** section ≈ 44 µm long, so that the pump (mid-gap) and the signal/idler (±1.25 THz) are all trapped in the defect. Tools: MPB 3D unit cell for the gap; coupled-mode theory (CMT) for the defect length; MEEP 2D (effective index n = 1.875) for the cavity modes, Q and fields. All in `defect_cavity_checks.py`.

| Check | Result |
|---|---|
| Band gap (MPB 3D) | converged < 1 % between resolution 32 and 48: **1.13 THz at 100 nm, 3.74 THz at 200 nm** corrugation. Mid-gap is at 1538 nm for a = 437 nm → period must grow ~0.8 % (a = 440.3 nm) to put it at 1547–1550 nm. Next band 16 THz above → single-mode. |
| Trapped modes | yes: four high-Q modes at 1540 / 1548 / 1557 / 1566 nm, spacing **1.0–1.3 THz**. Decay length into the grating 12 µm (pump) / 15 µm (signal). |
| F (mode-volume gain) | length-only estimate 13.5, but all three modes are **standing waves**, so the FWM overlap is ~4× smaller. **From the MEEP fields: F = 3.2–3.4 (Ey only), 2.4–2.5 (full vector) → F ≈ 3.** CMT gives 0.2–4.2 depending on sub-µm details of the defect length: the defect length is a sensitive design knob. |
| Radiation Q (MEEP 2D) | **abrupt** grating start: Q_rad 1×10⁴–7×10⁴. **20-period taper:** Q_rad 1.1×10⁵–2.1×10⁶, different for each mode (pump 1548 nm: 1.5–2.1×10⁶; 1540 nm: 1.1–1.5×10⁵; 1566 nm: 2.5–3.4×10⁵). Compare Q_int = 8×10⁵. |
| Frequency matching | f_s + f_i − 2f_p = **+36 GHz and −108 GHz** for the two tapered triplets (+39 GHz abrupt), i.e. 37–112 linewidths. CMT says a symmetric defect with the pump exactly at mid-gap can null it; not designed yet. |

**Net effect on the pump.** Each mode's total intrinsic loss grows by r = 1 + Q_int/Q_rad. At fixed η the required pump scales as κ_p·κ_s (the "powers scale as κ₀²" rule), so

  P_defect / P_ring ≈ r_p · √(r_s r_i) / F

| Design | P_defect / P_ring (if the mismatch were tuned to 0) |
|---|---|
| abrupt | **≈ 160–180** (much worse) |
| tapered, pump at 1557 nm | **1.1–3.3** |
| tapered, pump at 1548 nm | **1.5–2.9** |
| needed for 0.5× | Q_rad ≥ **4×10⁶ on all three modes** (r ≤ 1.2) + mismatch ≲ 1 GHz |

**Conclusion.** The defect does shrink the mode volume (F ≈ 3), but **as designed it costs pump**, because radiation at the defect adds loss to the signal/idler as well as the pump (it also lowers the escape efficiency of the squeezed light), and the frequencies are not matched. Today the MZI/two-point coupler (≈ 0.23× at 12 dB, no new loss mechanism) is the stronger lever. The defect route stays open only with a better taper (Q_rad > 4×10⁶ in **3D**, where vertical radiation is added) and a symmetric, frequency-matched defect; the two levers could in principle be combined.

Caveats: 2D effective-index model (no out-of-plane radiation → 3D Q will be lower); straight cavity (no ring curvature); Q values > 10⁶ are not well resolved in the simulated time; defect not tuned for matching.

### 6.2 Two-point (MZI-coupled) coupler in the engine (Nov 17)

**Model (PART 9b, `mzi_*()`).** The bus touches the ring at point 1 (rate κ₁), runs through a bypass (power transmission T_b = τ², phase ψ relative to the ring arc) and touches it again at point 2 (κ₂). For each mode (units of κ₀):

  Γ = (κ₀+κ₁+κ₂)/2 + τ√(κ₁κ₂)e^{iψ},  c_in = √κ₁ + τ√κ₂e^{iψ},  c_out = √κ₂ + τ√κ₁e^{iψ}

with a vacuum entering through the bypass loss. The squeezer is the same Bogoliubov/Langevin treatment as `ring_squeezer()`, but with separate Γ for signal and idler. Phases: ψ_p = π + ε, ψ_s = ε + δψ, ψ_i = ε − δψ (ε = heater error, δψ = bypass-delay error, ~3° per µm of bypass length). Design variables: √κ₁ = (u+v)/2, √κ₂ = (u−v)/2, where u² is the signal/idler coupling and v² the pump coupling (v = 1: pump critical). The optimiser chooses u, v (and which point is the strong one). **Validated:** κ₂ = 0 reproduces the plain ring of the engine exactly, and the ideal coupler reproduces `two_point_bound_dB()` exactly (4 digits).

**Fair comparison (pump / plain ring at its best η):**

| Case | 3 dB | 4 dB | 6 dB | 8 dB | 10 dB | 12 dB |
|---|---|---|---|---|---|---|
| Ideal (= the old bound) | 0.92× | 0.85× | 0.68× | 0.50× | 0.34× | 0.23× |
| Bypass loss 1% | 0.94× | 0.87× | 0.70× | 0.51× | 0.36× | 0.25× |
| Bypass loss 2% | 0.95× | 0.88× | 0.71× | 0.53× | 0.38× | 0.27× |
| Bypass loss 5% | 0.99× | 0.93× | 0.77× | 0.59× | 0.45× | 0.35× |
| Bypass loss 10% | 1.07× | 1.01× | 0.86× | 0.70× | 0.58× | 0.52× |
| Heater error 20° | 0.93× | 0.87× | 0.70× | 0.52× | 0.37× | 0.25× |
| Heater error 30° | 0.94× | 0.88× | 0.73× | 0.55× | 0.40× | 0.30× |
| Delay error 20° (~7 µm), heater re-tuned | 0.92× | 0.85× | 0.68× | 0.50× | 0.34× | 0.23× |
| Coupler dispersion 20% per 10 nm (η_s ≠ η_i) | 0.94× | 0.87× | 0.69× | 0.51× | 0.35× | 0.23× |
| **Realistic:** 1% bypass + 20%/10 nm + 10° delay + 10° heater | **0.95×** | **0.89×** | **0.72×** | **0.53×** | **0.38×** | **0.26×** |
| *PhC ring + pump reflector, ideal reflectors (§4.2)* | — | *1.62×* | — | *0.91×* | *0.61×* | *0.40×* |

At Vaidya's ring (η = 0.75, 100 mW), the realistic coupler gives 5.08 dB vs 5.23 dB ideal and 4.41 dB for the plain ring.

**What we learned:**
1. **Bypass loss is the main cost, and it is small.** A ~0.5 mm bypass at Q_int = 8×10⁵ (~0.46 dB/cm) loses ~0.5%. Even 5% still wins at every target. With bypass loss it is slightly better to put the **weak point first**, so that the signal leaking at the strong point does not cross the bypass.
2. **The heater error creates a coupling-induced resonance shift.** Im(Γ) = τ√(κ₁κ₂) sin ψ. At ψ ≠ 0, π the pump shifts one way and signal/idler the other, so the pair is no longer frequency matched (the earlier "20° → 1.2% pump leak" estimate ignored this). Still, 20° costs only ~0.02×, and retuning the laser hardly helps.
3. **A bypass-delay error shifts signal and idler in opposite directions.** This is cancelled exactly by choosing the detection sideband Ω (the model does this), so ±7 µm of length error costs nothing. At 13 µm, 12 dB becomes unreachable because η drops below 1 − 10^(−1.2) = 0.937.
4. **η_s ≠ η_i from coupler dispersion is harmless** up to ~20% per 10 nm.
5. **Correction to the earlier tuning plan** ("design the pump coupling slightly below the expected loss, raise it with the heater"): in the model the heater does not pay. Moving ψ_p away from π also lowers the signal coupling and adds frequency mismatch. If the real intrinsic loss is *lower* than designed (×0.8), nothing is lost (0.53× / 0.37× / 0.25× at 8 / 10 / 12 dB). If it is *higher* (×1.3), high targets become unreachable because η_s falls, which is true for any layout. So: **design for the pessimistic loss**.

**Not yet in the model:** Kerr (SPM/XPM), the spectrum Ω ≠ 0 beyond the sideband choice, the frequency dependence of ψ across a linewidth, and the real gap → κ₁, κ₂ numbers (needs a coupler cross-section calculation).

---

## 7. Other open items

| Item | Why it matters | Status |
|---|---|---|
| Two-point (MZI-coupled) coupler in the engine | It is the strongest lever | **done (§6.2)**; next: gap → κ₁, κ₂ for the real waveguide, and Kerr |
| Frequency matching at the split pump mode | The laser sits ~4.7 linewidths from the unperturbed resonance; the model assumes signal + idler are matched to it (NIST uses PhC splitting exactly to fix frequency mismatch) | not modelled |
| Kerr (SPM/XPM) and the spectrum Ω ≠ 0 | Winning points are near threshold, 0.4–2 W at Q_int = 8×10⁵ | not modelled |
| Source of the x formula (`x_plain`) | Reproduces Vaidya's ~4 dB at 100 mW, but not checked against their supplement | to verify |

**Questions for Eran:** which squeezing target we design for (PhC + pump reflector only pays above ~8–9 dB); what PhC he has in mind (uniform corrugation or a localized/defect mode) and whether signal and idler should also be localized; realistic Bragg reflector loss (< 0.5%?); whether to model the MZI coupler in parallel.

**Prior art note:** the two-point coupler is a microring coupled through a Mach–Zehnder interferometer, a known structure (Barbarossa 1995; Tison et al. 2017; Vernon et al. 2017) used for photon-pair sources to set pump and signal coupling separately. No CW two-mode squeezing use was found. PhC rings are used for squeezing by DESY (Ulanov et al. 2025, ~7.8 dB on-chip, corrugation to suppress parasitic processes) and for OPO frequency matching by NIST (Lu et al. 2022). Using the split pump mode + a bus reflector as a pump coupling knob was not found in the literature.

---

## 8. How to run MPB

MPB has no GUI; it is a Python library (`from meep import mpb`). It is not available natively on Windows (use WSL2 with Ubuntu). On Linux/macOS with Miniforge:

```
conda create -n mp -c conda-forge python=3.11 pymeep numpy matplotlib
conda activate mp
python xsec.py
```
