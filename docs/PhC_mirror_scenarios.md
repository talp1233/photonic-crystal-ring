# PhC ring + mirrors: three layouts, fair comparison

*Prepared for the meeting with Eran, Nov 2026 · Tal. Built on Yaniv's `ring_three_configs.py` engine (unchanged); script: `phc_mirror_scenarios.py`.*

## The message

A PhC ring with a pump mirror on the right is the first layout that **beats the plain ring** in the fair comparison, but only at high targets: above about 7.5 dB with ideal mirrors, or about 9 dB with 2% mirrors. At 10 dB it needs 0.61× the pump (ideal) or 0.87× (2%). At 12 dB the figures are 0.40× and 0.75×.

**Why it works:** the mirror phase sets the pump's effective coupling independently of the coupling gap. So the pump can be critically coupled (95–100% of it used) while signal and idler stay overcoupled through the gap (best η ≈ 0.9). This is the split "pump critical / squeezing overcoupled" that Eran asked for.

**What it costs:** the pump is a standing wave, so its energy splits equally between CW and CCW. That gives two half-strength squeezers on two ports. As a result the layout loses below about 8 dB (1.3–1.8×), and it stays about 1.7–1.8× above the ideal bound, where all the pump drives one direction (e.g. the two-point coupler).

**Sending one half back through the ring (B, C) adds essentially nothing.** At the operating points the crossing squeezer has t² ≤ 0.03, so the returning squeezed light is erased. C also adds a second mirror loss.

## The three layouts

Pump enters from the left in all three. The PhC splits **only the pump resonance**, by 10 Vaidya linewidths (≈10 GHz, fixed in absolute terms). Signal and idler remain travelling waves. Mirrors use Yaniv's placeholder: R = 97%, loss 2%, T = 1% in the reflect band, and 98% transmission in the pass band.

| Layout | Mirrors | Useful port |
|---|---|---|
| **A** | pump mirror, right | either port (they are nearly equal; the better one is used) |
| **B** | full mirror (pump + squeezing), right | left |
| **C** | pump mirror right + squeezing mirror left | right |

For each layout and η, the model chooses:
- the laser detuning (to sit on a split mode)
- the pump-mirror phase (to maximize the pump used)
- the phase of any returning squeezed light: best case. For B we also show "in-step", meaning returning pair phase = pump phase. This is the default for a non-dispersive path, since the mirror distance cannot tune it.

## 1. At Vaidya's ring (η = 0.75, 100 mW, 2% mirrors)

| Layout | Out left (dB) | Out right (dB) | Pump used (mW) | Pump out left / right (mW) | Lost in mirrors (mW) | x CCW / CW |
|---|---|---|---|---|---|---|
| 1. Ring alone | 0 | 4.41 | 75.0 | 0 / 25.0 | 0 | 0.44 / 0 |
| 2. Ring + pump mirror | 1.32 | 4.26 | 93.2 | 6.1 / 0.3 | 0.5 | 0.44 / 0.11 |
| 3. PhC ring alone | 1.36 | 1.37 | 37.6 | 56.1 / 6.3 | 0 | 0.11 / 0.11 |
| **A** PhC + pump mirror | 3.28 | 2.89 | 95.4 | 0.1 / 1.5 | 2.9 | 0.26 / 0.30 |
| **B** PhC + full mirror | 3.39 (in-step 3.24, worst 2.98) | 0.02 | 95.4 | 0.1 / 1.5 | 2.9 | 0.26 / 0.30 |
| **C** PhC + pump + sq. mirrors | 0.02 | 3.00 (worst 2.36) | 93.5 | 0.1 / 1.4 | 4.9 | 0.26 / 0.29 |

With ideal mirrors, A, B and C use 100% of the pump, and give 3.4 / 3.5 / 3.2 dB at the useful port. For reference, the ideal bound at this point (all pump used, one direction) gives 5.2 dB.

**The point:** the PhC with a mirror fixes the pump (95% used and almost nothing back to the laser), but at η = 0.75 each squeezer has only about 0.3 drive, so it still loses to the plain ring.

## 2. Same power (100 mW), different η, 2% mirrors

All values are squeezing at the useful port, in dB.

| η | 1. Ring | 3. PhC | **A** left / right | **B** | **C** | Bound (no mirrors) |
|---|---|---|---|---|---|---|
| 0.60 | 3.96 | 1.95 | 3.23 / 3.02 | 3.39 | 3.10 | 3.97 |
| 0.70 | 4.62 | 1.65 | 3.41 / 3.07 | 3.41 | 3.04 | 4.94 |
| 0.75 | 4.41 | 1.37 | 3.28 / 2.89 | 3.39 | 3.00 | 5.23 |
| 0.80 | 3.66 | 1.03 | 2.66 / 2.89 | 3.08 | 3.16 | 5.18 |
| 0.85 | 2.52 | 0.67 | 2.13 / 2.40 | 2.87 | 2.92 | 4.63 |
| 0.90 | 1.30 | 0.35 | 1.45 / 1.72 | 2.36 | 2.36 | 3.52 |
| 0.95 | 0.37 | 0.11 | 0.67 / 0.87 | 1.33 | 1.29 | 1.93 |

At fixed η, the PhC layouts overtake the plain ring above η ≈ 0.85. The plain ring loses fast at high η because most of its pump passes by.

## 3. Same power, each at its best η

Each cell is dB at the useful port, with the best η in brackets.

| Pump (mW) | Plain ring | A ideal | B ideal | C ideal | A 2% | B 2% | C 2% |
|---|---|---|---|---|---|---|---|
| 50 | 3.17 (0.63) | 2.07 (0.61) | 2.08 (0.69) | 2.08 (0.67) | 2.03 (0.60) | 2.03 (0.60) | 1.95 (0.62) |
| 100 | 4.62 (0.70) | 3.50 (0.70) | 3.51 (0.75) | 3.49 (0.74) | 3.41 (0.69) | 3.41 (0.68) | 3.26 (0.66) |
| 200 | 6.16 (0.78) | 5.44 (0.77) | 5.47 (0.82) | 5.44 (0.81) | 5.26 (0.75) | 5.29 (0.75) | 5.03 (0.76) |
| 400 | 7.72 (0.84) | 7.89 (0.86) | 7.89 (0.86) | 7.85 (0.88) | 7.17 (0.85) | 7.14 (0.79) | 7.11 (0.83) |
| 800 | 9.22 (0.89) | 10.67 (0.92) | 10.67 (0.92) | 10.59 (0.93) | 9.34 (0.91) | 9.34 (0.90) | 9.29 (0.89) |

## 4. The fair comparison: pump needed for a target, relative to the plain ring

Intrinsic Q is fixed at 8×10⁵, and every layout uses its best η. Values below 1× mean the layout needs less pump than the plain ring.

| Target | Plain ring: mW (best η) | Bound (two-point coupler, ideal) |
|---|---|---|
| 3 dB | 46 (0.62) | 0.92× |
| 4 dB | 75 (0.67) | 0.85× |
| 6 dB | 186 (0.77) | 0.68× |
| 8 dB | 454 (0.85) | 0.50× |
| 10 dB | 1113 (0.90) | 0.34× |
| 12 dB | 2748 (0.94) | 0.23× |

**Mirror loss 0% (ideal)**

| Target | 4. Ring + 2 mirrors | 5. Ring + full | **A** | **B** (best / in-step) | **C** |
|---|---|---|---|---|---|
| 3 dB | 1.00× | 3.29× | 1.76× (η 0.67) | 1.75× / 1.76× | 1.76× |
| 4 dB | 0.99× | 3.13× | 1.62× (η 0.72) | 1.61× / 1.62× | 1.62× |
| 6 dB | 0.95× | 2.64× | 1.27× (η 0.80) | 1.26× / 1.27× | 1.27× |
| 8 dB | 0.88× | 2.06× | 0.91× (η 0.86) | 0.91× / 0.91× | 0.92× |
| 10 dB | 0.88× | 1.67× | **0.61×** (η 0.91) | 0.61× / 0.61× | 0.62× |
| 12 dB | 0.92× | 1.43× | **0.40×** (η 0.95) | 0.40× / 0.40× | 0.41× |

**Mirror loss 0.5%**

| Target | 4. Ring + 2 mirrors | 5. Ring + full | **A** | **B** (best / in-step) | **C** |
|---|---|---|---|---|---|
| 3 dB | 1.01× | 3.32× | 1.78× | 1.77× / 1.78× | 1.80× |
| 4 dB | 1.01× | 3.16× | 1.65× | 1.63× / 1.66× | 1.67× |
| 6 dB | 0.98× | 2.67× | 1.29× | 1.29× / 1.29× | 1.33× |
| 8 dB | 0.93× | 2.08× | 0.93× | 0.93× / 0.93× | 0.97× |
| 10 dB | 0.98× | 1.69× | 0.64× | 0.64× / 0.64× | 0.68× |
| 12 dB | 1.09× | 1.44× | 0.43× | 0.42× / 0.42× | 0.47× |

**Mirror loss 2% (placeholder)**

| Target | 4. Ring + 2 mirrors | 5. Ring + full | **A** | **B** (best / in-step) | **C** |
|---|---|---|---|---|---|
| 3 dB | 1.07× | 3.42× | 1.82× | 1.82× / 1.83× | 1.93× |
| 4 dB | 1.08× | 3.25× | 1.69× | 1.69× / 1.70× | 1.81× |
| 6 dB | 1.09× | 2.74× | 1.47× | 1.35× / 1.53× | 1.50× |
| 8 dB | 1.13× | 2.14× | 1.15× | 1.15× / 1.17× | 1.17× |
| 10 dB | 1.38× | 1.74× | 0.87× | 0.87× / 0.87× | 0.91× |
| 12 dB | 2.00× | 1.50× | 0.75× | 0.67× / 0.67× | 0.76× |

**Check:** the columns for scenarios 4 and 5 reproduce Yaniv's eta_sweep mirror-loss table to within about 0.03.

## Conclusions

1. **PhC + pump mirror (A) is the first layout that wins, and it wins at high targets.** The pump-mirror phase acts as an independent pump coupling, which gives critical coupling for the pump with overcoupled squeezing.
2. **The remaining gap to the bound (about 1.7–1.8×) is the standing wave.** Each direction gets half the pump energy. To close it, all the pump must drive one direction: for example, also split signal and idler (open question), or use the two-point coupler.
3. **B and C do not improve on A.** Returning squeezed light is erased when it crosses the ring, since t² ≈ 0 at the operating points. A is the simplest layout and needs only one mirror.
4. **Mirror loss still matters, but less than for the plain-ring layouts.** In A, the mirror loss sits inside the pump's composite cavity (extra pump loss) and on the right port. Going from 2% to 0.5% changes 10 dB from 0.87× to 0.64×.

## Caveats

- **Frequency matching.** The laser sits on a split pump mode, about 4.7 linewidths below the unperturbed resonance. The model assumes signal + idler are frequency-matched to that shifted pump. NIST uses exactly this splitting to correct frequency mismatch in PhC-ring OPOs (Lu et al., Opt. Lett. 2022), so the ring dispersion must be designed for it.
- **Splitting.** It is fixed at about 10 GHz absolute, so at η = 0.9 it is only 4 loaded linewidths. Sensitivity for A with ideal mirrors:

  | Splitting | 10 dB | 12 dB |
  |---|---|---|
  | 0.5× | 600 mW | 2607 mW |
  | 1× (as used) | 682 mW | 1248 mW |
  | 2× | 722 mW | 1173 mW |

- **High power, near threshold, no Kerr.** The winning targets need 0.4–2 W, and the best points sit near threshold. Kerr and thermal effects are neglected.
- **Mirror distance** is assumed short: no extra delay, single-frequency phase.
- **Ω = 0 only.** Only the centre of the squeezing spectrum is modelled.

## Questions for Eran

1. Which resonances should the PhC split: pump only, or signal and idler too? And how large a splitting is realistic?
2. Is the frequency matching at the split pump mode acceptable as a design constraint, i.e. dispersion engineered so that the shifted pump is matched?
3. Which squeezing target are we designing for? A only pays off above about 8–9 dB.
4. Realistic loss of a waveguide Bragg pump reflector: below 0.5%?
5. Given that A ≈ B ≈ C, should we drop B and C, and put the next effort into a PhC that splits all three resonances, or into the two-point coupler?
