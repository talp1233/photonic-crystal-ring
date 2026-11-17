#!/usr/bin/env python3
# ============================================================
#  Level 2 ground rules + three configurations
# ============================================================
#
#  Layout (left to right on the bus waveguide):
#
#      pump in  -->   [left mirror?]   ( ring )   [right mirror?]   -->
#
#  Directions inside the ring:
#      CCW mode = driven by light moving RIGHT on the bus (the input pump).
#                 Its squeezed pairs leave the ring moving RIGHT.
#      CW  mode = driven by light moving LEFT on the bus (anything a right
#                 mirror sends back). Its squeezed pairs leave moving LEFT.
#
#  Ring types:
#      regular ring : CW and CCW never talk to each other.
#      PhC ring     : a grating in the ring couples CW <-> CCW at the PUMP
#                     mode only (mode splitting). Signal/idler modes stay
#                     ordinary travelling waves.
#
#  The three simulations:
#      Sim 1: PhC ring vs regular ring, no mirrors.
#      Sim 2: regular ring + pump mirror on the RIGHT + squeezing mirror on
#             the LEFT  -> pump ends up on the left, squeezing on the right.
#      Sim 3: regular ring + full mirror on the RIGHT (pump + squeezing)
#             -> everything ends up on the left.
#
#  Mirrors have LOSS: at every band  R + T + A = 1.
#
#  Same published ring as before (Xanadu, Vaidya 2020). Closed forms:
#      pump      : coupled-mode theory, steady state (solved exactly)
#      squeezing : Collett-Gardiner cavity input-output at the centre of
#                  the squeezing spectrum (Omega = 0), written as a
#                  Bogoliubov map so squeezed light can be fed INTO a
#                  second cavity (needed for configurations 2 and 3).
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent

# ------------------------------------------------------------
# STEP 0 -- fixed numbers
# ------------------------------------------------------------
# Ring (Vaidya 2020) -- fabrication-set, never swept
Q          = 2.0e5     # loaded Q (same for pump and squeezed modes)
eta        = 0.75      # escape efficiency (pump and squeezed modes)
wavelength = 1550.2e-9
hbar, c    = 1.054e-34, 3.0e8
ring_radius, FSR = 120e-6, 190e9

w      = 2*np.pi*c/wavelength
length = 2*np.pi*ring_radius
n_g    = c/(FSR*length)
Lambda = hbar*w*(c/n_g)**2*1.0/length      # nonlinear strength (gamma_NL = 1)
half_linewidth = w/(2*Q)

# Mirrors -- PLACEHOLDERS until we have a real grating design.
# "reflect band": the band the mirror is built to bounce back
# "pass band"   : the bands it is supposed to let through
MIRROR_R      = 0.97   # power reflectivity in its reflect band
MIRROR_LOSS   = 0.02   # absorption + scattering, any band (0.09 dB)
# (so in the reflect band T = 1 - 0.97 - 0.02 = 0.01,
#  and in the pass band T = 1 - 0.02 = 0.98)

# PhC ring: CW<->CCW coupling at the pump mode, in units of the linewidth.
# beta = 5 -> splitting 2*beta = 10 linewidths (~10 GHz), "strong" splitting.
PHC_BETA = 5.0


# ------------------------------------------------------------
# STEP 1 -- the drive strength from pump power (unchanged from before)
# ------------------------------------------------------------
def x_plain(power_mW):
    """Drive strength x of a REGULAR ring with all pump in one direction.
    Same formula as ring_vs_reflector.py (x = 1 is threshold)."""
    pump_inside = 4*eta*Q*(power_mW/1000)/(hbar*w**2)
    return Lambda*pump_inside/half_linewidth


# ------------------------------------------------------------
# STEP 2 -- classical pump: where does the pump go?
# ------------------------------------------------------------
# Units: linewidth kappa = 1, input amplitude s = 1 (so |s|^2 = P_in).
# Unknowns: a1 = CCW pump amplitude, a2 = CW pump amplitude,
#           SR = field hitting the ring from the left (moving right),
#           SL = field hitting the ring from the right (moving left).
#
#   (1/2 - i*D) a1 - i*beta a2 = sqrt(eta) SR       ring, CCW
#   (1/2 - i*D) a2 - i*beta a1 = sqrt(eta) SL       ring, CW
#   OR = SR - sqrt(eta) a1     light leaving the ring to the right
#   OL = SL - sqrt(eta) a2     light leaving the ring to the left
#   SR = tL * s                (left mirror never reflects pump here)
#   SL = rR * exp(i*phiR) * OR (right mirror sends some back)
#
# A regular ring is just beta = 0.
def solve_pump(beta=0.0, detuning=0.0, left=None, right=None, phiR=0.0):
    TL = 1.0 if left is None else left["T_pump"]
    AL = 0.0 if left is None else left["A_pump"]
    RR = 0.0 if right is None else right["R_pump"]
    TR = 1.0 if right is None else right["T_pump"]
    AR = 0.0 if right is None else right["A_pump"]
    se, D, rRe = np.sqrt(eta), detuning, np.sqrt(RR)*np.exp(1j*phiR)
    tL = np.sqrt(TL)

    # unknowns [a1, a2, SR, SL]
    M = np.array([[0.5-1j*D, -1j*beta, -se, 0],
                  [-1j*beta, 0.5-1j*D, 0, -se],
                  [0, 0, 1, 0],
                  [rRe*se, 0, -rRe, 1]], dtype=complex)
    b = np.array([0, 0, tL, 0], dtype=complex)
    a1, a2, SR, SL = np.linalg.solve(M, b)
    OR, OL = SR - se*a1, SL - se*a2

    out = dict(a_ccw=a1, a_cw=a2,
               N_ccw=abs(a1)**2, N_cw=abs(a2)**2,
               pump_right=TR*abs(OR)**2,              # exits on the right
               pump_left=TL*abs(OL)**2,               # exits on the left
               used_in_ring=(1-eta)*(abs(a1)**2+abs(a2)**2),
               lost_in_mirrors=AL*(1+abs(OL)**2) + AR*abs(OR)**2)
    out["energy_check"] = (out["pump_right"] + out["pump_left"]
                           + out["used_in_ring"] + out["lost_in_mirrors"])
    return out


def best_pump(beta, left=None, right=None, scan_phase=False):
    """Pick the laser detuning (and right-mirror phase, if asked) that puts
    the most pump to use in the ring. For beta = 0 that is simply D = 0."""
    dets = [0.0] if beta == 0 else np.linspace(-1.5*beta, 1.5*beta, 601)
    phis = np.linspace(0, 2*np.pi, 181) if scan_phase else [0.0]
    best = None
    for D in dets:
        for p in phis:
            s = solve_pump(beta, D, left, right, p)
            if best is None or s["used_in_ring"] > best["used_in_ring"]:
                best, best["D"], best["phiR"] = s, D, p
    return best


# ------------------------------------------------------------
# STEP 3 -- squeezed light as a Bogoliubov map
# ------------------------------------------------------------
# Every field is written as   f = sum_k alpha_k c_k + beta_k c_k^dagger
# where the c_k are independent vacuum inputs (bus, ring loss, mirror loss).
class Vacua:
    def __init__(self, n_max=60):
        self.n, self.n_max = 0, n_max
    def new(self):
        al, be = np.zeros(self.n_max, complex), np.zeros(self.n_max, complex)
        al[self.n] = 1.0
        self.n += 1
        return (al, be)

def dag(f):            return (np.conj(f[1]), np.conj(f[0]))
def scale(f, z):       return (z*f[0], z*f[1])
def add(*fs):          return (sum(f[0] for f in fs), sum(f[1] for f in fs))

def through_loss(f, keep, vac):
    """Pass a field through something that keeps a fraction 'keep' of the
    power; the rest is replaced by fresh vacuum."""
    return add(scale(f, np.sqrt(keep)), scale(vac.new(), np.sqrt(1-keep)))

def ring_squeezer(fs, fi, x, theta, vac):
    """One direction of the ring at line centre (Collett-Gardiner).
    fs, fi : signal and idler arriving on the bus
    x      : drive strength (x = 1 is threshold), theta : pump phase
    Returns the signal and idler leaving the ring on the bus."""
    G = x/2
    A = np.array([[0.5, -G*np.exp(1j*theta)],
                  [-G*np.exp(-1j*theta), 0.5]])
    Ainv = np.linalg.inv(A)
    ls, li = vac.new(), vac.new()                         # ring-loss vacua
    drive = [add(scale(fs, np.sqrt(eta)), scale(ls, np.sqrt(1-eta))),
             add(scale(dag(fi), np.sqrt(eta)), scale(dag(li), np.sqrt(1-eta)))]
    cav = [add(scale(drive[0], Ainv[r, 0]), scale(drive[1], Ainv[r, 1]))
           for r in range(2)]
    out_s  = add(fs,      scale(cav[0], -np.sqrt(eta)))
    out_id = add(dag(fi), scale(cav[1], -np.sqrt(eta)))
    return out_s, dag(out_id)

def squeezing_dB(bs, bi):
    """Best two-mode squeezing (dB, positive = below shot noise) and the
    anti-squeezing, from the covariance matrix of the signal/idler pair."""
    rows = []
    for b in (bs, bi):
        u, v = (b[0]+b[1])/2, 1j*(b[0]-b[1])/2
        rows.append(np.concatenate([2*u.real, 2*v.real]))   # q quadrature
        rows.append(np.concatenate([2*u.imag, 2*v.imag]))   # p quadrature
    R = np.array(rows)
    ev = np.linalg.eigvalsh(R @ R.T)
    return -10*np.log10(ev[0]), 10*np.log10(ev[-1])

def made_in_ring_dB(x):
    """Squeezing the process itself makes (escape efficiency set to 1)."""
    return -20*np.log10((1-x)/(1+x)) if x < 1 else np.inf


# ------------------------------------------------------------
# STEP 4 -- one function that runs any configuration
# ------------------------------------------------------------
def make_mirror(reflects_pump, reflects_squeezing):
    m = {}
    for band, refl in (("pump", reflects_pump), ("sq", reflects_squeezing)):
        R = MIRROR_R if refl else 0.0
        m["R_"+band], m["A_"+band] = R, MIRROR_LOSS
        m["T_"+band] = 1 - R - MIRROR_LOSS
    return m

def run(power_mW, ring="regular", left=None, right=None,
        sq_phase=0.0, pump=None):
    """left/right : mirror dicts or None.
    sq_phase    : extra round-trip phase of the squeezed pair (signal+idler)
                  relative to what the pump picks up. 0 = the two stay in step.
    Returns pump bookkeeping + squeezing on each side."""
    beta = PHC_BETA if ring == "PhC" else 0.0
    if pump is None:
        pump = best_pump(beta, left, right,
                         scan_phase=(right is not None and beta > 0))
    x0 = x_plain(power_mW)
    N_ref = 4*eta                              # plain ring, all pump one way
    x_ccw = x0*pump["N_ccw"]/N_ref
    x_cw  = x0*pump["N_cw"]/N_ref
    th_ccw, th_cw = 2*np.angle(pump["a_ccw"]), 2*np.angle(pump["a_cw"])

    vac = Vacua()
    # --- squeezed light that ENTERS the ring from the left, moving right
    fs, fi = vac.new(), vac.new()
    if left is not None:            # it came through the left mirror
        fs, fi = through_loss(fs, left["T_sq"], vac), through_loss(fi, left["T_sq"], vac)

    # Order of events depends on which mirror bounces squeezed light.
    if left is not None and left["R_sq"] > 0:
        # CONFIG 2: the CW pairs (made by the recycled pump) go left, bounce
        # off the squeezing mirror, and come back into the CCW mode.
        gs, gi = vac.new(), vac.new()                     # vacuum from the right
        cw_s, cw_i = ring_squeezer(gs, gi, x_cw, th_cw, vac)
        # leftward CW output hits the left mirror
        back_s = through_loss(cw_s, left["R_sq"], vac)
        back_i = through_loss(cw_i, left["R_sq"], vac)
        ph = np.exp(1j*(sq_phase/2))
        back_s, back_i = scale(back_s, ph), scale(back_i, ph)
        # what arrives moving right = reflected CW light + leak of the outside vacuum
        fs = add(back_s, scale(vac.new(), np.sqrt(left["T_sq"])))
        fi = add(back_i, scale(vac.new(), np.sqrt(left["T_sq"])))
        ccw_s, ccw_i = ring_squeezer(fs, fi, x_ccw, th_ccw, vac)
        right_s, right_i = ccw_s, ccw_i
        if right is not None:       # pump mirror on the right: pairs pass with loss
            right_s = through_loss(ccw_s, right["T_sq"], vac)
            right_i = through_loss(ccw_i, right["T_sq"], vac)
        left_s  = through_loss(cw_s, left["T_sq"], vac)
        left_i  = through_loss(cw_i, left["T_sq"], vac)
    else:
        ccw_s, ccw_i = ring_squeezer(fs, fi, x_ccw, th_ccw, vac)
        if right is not None and right["R_sq"] > 0:
            # CONFIG 3: CCW pairs go right, bounce, come back into the CW mode.
            # The pump took round-trip phase phiR; the pair takes 2*phiR - sq_phase.
            ph = np.exp(1j*(pump["phiR"] - sq_phase/2))
            back_s = scale(through_loss(ccw_s, right["R_sq"], vac), ph)
            back_i = scale(through_loss(ccw_i, right["R_sq"], vac), ph)
            right_s = through_loss(ccw_s, right["T_sq"], vac)
            right_i = through_loss(ccw_i, right["T_sq"], vac)
        else:
            back_s, back_i = vac.new(), vac.new()
            if right is not None:   # pump-only mirror: pairs just pass (with loss)
                right_s = through_loss(ccw_s, right["T_sq"], vac)
                right_i = through_loss(ccw_i, right["T_sq"], vac)
            else:
                right_s, right_i = ccw_s, ccw_i
        cw_s, cw_i = ring_squeezer(back_s, back_i, x_cw, th_cw, vac)
        left_s, left_i = cw_s, cw_i
        if left is not None:
            left_s = through_loss(left_s, left["T_sq"], vac)
            left_i = through_loss(left_i, left["T_sq"], vac)

    sqR, antiR = squeezing_dB(right_s, right_i)
    sqL, antiL = squeezing_dB(left_s, left_i)
    return dict(pump=pump, x_ccw=x_ccw, x_cw=x_cw,
                made_ccw=made_in_ring_dB(x_ccw), made_cw=made_in_ring_dB(x_cw),
                sq_right=sqR, anti_right=antiR, sq_left=sqL, anti_left=antiL)


# ============================================================
# RUN
# ============================================================
def set_eta(value):
    """Only used for the 'what if the ring were more over-coupled' parts."""
    global eta
    eta = value

PUMP_MIRROR = make_mirror(reflects_pump=True,  reflects_squeezing=False)
SQ_MIRROR   = make_mirror(reflects_pump=False, reflects_squeezing=True)
FULL_MIRROR = make_mirror(reflects_pump=True,  reflects_squeezing=True)

# layout name -> (ring, mirrors, port where the useful squeezing comes out)
LAYOUTS = {
    "alone":  ("regular", {},                                       "right"),
    "pumpR":  ("regular", dict(right=PUMP_MIRROR),                  "right"),  # ground rules
    "sim1":   ("PhC",     {},                                       "right"),
    "sim2":   ("regular", dict(left=SQ_MIRROR, right=PUMP_MIRROR),  "right"),
    "sim3":   ("regular", dict(right=FULL_MIRROR),                  "left"),
    "aside":  ("PhC",     dict(right=FULL_MIRROR),                  "left"),   # not requested
}

def solve_layouts():
    """The pump is linear, so each layout's pump is solved once (any power)."""
    out = {}
    for name, (ring, mir, _) in LAYOUTS.items():
        beta = PHC_BETA if ring == "PhC" else 0.0
        out[name] = best_pump(beta, scan_phase=(beta > 0 and "right" in mir), **mir)
    return out

def R(P, name, pumps, ph=0.0):
    ring, mir, _ = LAYOUTS[name]
    return run(P, ring, sq_phase=ph, pump=pumps[name], **mir)

def phase_scan(P, name, pumps, n=181):
    side = LAYOUTS[name][2]
    ph = np.linspace(0, 2*np.pi, n)
    return ph, np.array([R(P, name, pumps, q)["sq_"+side] for q in ph])

def best_phase(P, name, pumps):
    """Squeezed light that returns to the ring has a free phase -> pick the best."""
    if not any(LAYOUTS[name][1].get(k) is not None and LAYOUTS[name][1][k]["R_sq"] > 0
               for k in ("left", "right")):
        return 0.0, None
    ph, v = phase_scan(P, name, pumps)
    return ph[np.argmax(v)], (v.min(), v.max())

def power_for(target_dB, name, pumps, ph):
    side = LAYOUTS[name][2]
    for P in np.arange(0.5, 400, 0.5):
        r = R(P, name, pumps, ph)
        if max(r["x_ccw"], r["x_cw"]) >= 1:
            return np.nan
        if r["sq_"+side] >= target_dB:
            return P
    return np.nan


if __name__ == "__main__":
    P0 = 100.0                                   # mW, Vaidya operating point
    pumps = solve_layouts()

    # ---------------- 1. validity checks ----------------
    print("="*76)
    print("1. CHECKS")
    print("="*76)
    x = x_plain(P0)
    cg = -10*np.log10(1 - eta*4*x/(1+x)**2)
    print(f"drive strength x at {P0:.0f} mW, regular ring : {x:.3f}  (threshold ~{P0/x:.0f} mW)")
    print(f"engine vs Collett-Gardiner closed form    : {R(P0,'alone',pumps)['sq_right']:.4f}"
          f" vs {cg:.4f} dB")
    print(f"no pump -> vacuum (expect 0 dB)           : {abs(R(1e-9,'alone',pumps)['sq_right']):.1e} dB")
    print(f"escape-limited ceiling                    : {-10*np.log10(1-eta):.1f} dB  (paper: ~75% escape)")
    print(f"out at {P0:.0f} mW, ring alone                : {R(P0,'alone',pumps)['sq_right']:.1f} dB"
          "  (paper: ~4 dB on-chip)")
    for name, pm in pumps.items():
        print(f"pump energy budget closes, {name:6s}         : {pm['energy_check']:.6f}")
    print()

    # ---------------- 2. ground rules ----------------
    print("="*76)
    print("2. GROUND RULES  (regular ring + pump mirror on the right)")
    print("="*76)
    used = [solve_pump(0, 0, None, PUMP_MIRROR, ph)["used_in_ring"]
            for ph in np.linspace(0, 2*np.pi, 13)]
    print(f"pump put to use vs mirror distance/phase  : min {min(used):.4f}, max {max(used):.4f}"
          "  -> distance does not matter (no cavity)")
    a, b = R(P0, "alone", pumps), R(P0, "pumpR", pumps)
    xb = x*(1 + 0.99*(1-2*eta)**2)
    old = -10*np.log10(1 - eta*4*xb/(1+xb)**2)
    print(f"pump put to use: ring alone {a['pump']['used_in_ring']:.3f} -> with mirror "
          f"{b['pump']['used_in_ring']:.3f}  (CW and CCW add as powers)")
    print(f"drive strength:  CCW {b['x_ccw']:.3f} (unchanged), CW {b['x_cw']:.3f} (new)")
    print(f"squeezing RIGHT : {a['sq_right']:.2f} dB alone -> {b['sq_right']:.2f} dB with mirror"
          "   (only the mirror's pass-band loss)")
    print(f"squeezing LEFT  : {a['sq_left']:.2f} dB alone -> {b['sq_left']:.2f} dB with mirror"
          "   (new, separate squeezer)")
    print(f"old 'boost' model would have said         : {old:.2f} dB on the right  <-- not physical")
    print()

    # ---------------- 3. the three simulations at P0 ----------------
    rows_ = [("alone", "reference: regular ring alone"),
             ("sim1",  "Sim 1: PhC ring, no mirrors"),
             ("sim2",  "Sim 2: regular, pump mirror R + squeeze mirror L"),
             ("sim3",  "Sim 3: regular, full mirror R"),
             ("aside", "aside: PhC ring + full mirror R")]
    best_ph = {}
    print("="*76)
    print(f"3. THREE SIMULATIONS at {P0:.0f} mW   (dB below shot noise; IN = made in ring)")
    print("="*76)
    for name, title in rows_:
        ph, rng = best_phase(P0, name, pumps)
        best_ph[name] = ph
        r = R(P0, name, pumps, ph)
        p = r["pump"]
        Pn = power_for(3.0, name, pumps, ph)
        print(title)
        print(f"   drive x  CCW {r['x_ccw']:.3f}  CW {r['x_cw']:.3f}   | made in ring  CCW {r['made_ccw']:.2f}"
              f"  CW {r['made_cw']:.2f} dB")
        print(f"   OUT right {r['sq_right']:.2f} dB   OUT left {r['sq_left']:.2f} dB"
              + (f"   (free phase: {rng[0]:.2f} .. {rng[1]:.2f} dB, best shown)" if rng else ""))
        print(f"   pump: used {p['used_in_ring']:.0%}  exits left {p['pump_left']:.0%}  "
              f"exits right {p['pump_right']:.0%}  lost in mirrors {p['lost_in_mirrors']:.0%}")
        print(f"   pump for 3 dB at the {LAYOUTS[name][2]} port: "
              + (f"{Pn:.1f} mW" if np.isfinite(Pn) else "not reached below threshold"))
    print()

    # ---------------- 4. what if the ring were more over-coupled? ----------------
    print("="*76)
    print("4. WHAT IF: same loaded Q, higher escape efficiency  (pump for 3 dB, best phase)")
    print("="*76)
    for e in (0.75, 0.90, 0.95):
        set_eta(e)
        pe = solve_layouts()
        need = {n: power_for(3.0, n, pe, best_phase(P0, n, pe)[0])
                for n in ("alone", "sim2", "sim3", "aside")}
        print(f"eta = {e:.2f}:  ring alone {need['alone']:6.1f} | sim 2 {need['sim2']:6.1f} | "
              f"sim 3 {need['sim3']:6.1f} | aside PhC+full {need['aside']:6.1f}  mW")
    set_eta(0.75)
    print()

    # ---------------- 5. side-by-side comparison of the five scenarios ----------------
    FIVE = [("alone", "1. Regular ring"),
            ("pumpR", "2. Ring + pump Bragg (R)"),
            ("sim1",  "3. PhC ring, no reflectors"),
            ("sim2",  "4. Ring + pump Bragg (R) + squeeze Bragg (L)"),
            ("sim3",  "5. Ring + full mirror (R)")]
    print("="*104)
    print(f"5. COMPARISON AT {P0:.0f} mW INPUT   (squeezing in dB below shot noise; "
          "'in ring' = made by each direction)")
    print("="*104)
    print(f"{'scenario':46s}| {'LEFT':>5} | {'in ring CCW':>11} {'in ring CW':>10} | {'RIGHT':>5} |"
          f" {'pump used':>9} {'out L':>6} {'out R':>6} {'mirrors':>7}")
    print("-"*104)
    summary = []
    for name, title in FIVE:
        ph, _ = best_phase(P0, name, pumps)
        r = R(P0, name, pumps, ph)
        p = r["pump"]
        summary.append((title, r))
        print(f"{title:46s}| {r['sq_left']:5.2f} | {r['made_ccw']:11.2f} {r['made_cw']:10.2f} | "
              f"{r['sq_right']:5.2f} | {P0*p['used_in_ring']:6.1f} mW {P0*p['pump_left']:5.1f}"
              f" {P0*p['pump_right']:6.1f} {P0*p['lost_in_mirrors']:6.1f}")
    print("pump columns in mW; scenarios 4-5 at the best phase of the returning squeezed light\n")

    # ============================================================
    # FIGURES
    # ============================================================
    BLUE, ORANGE = "#2a78d6", "#eb6834"      # regular, PhC
    INK, MUTED = "#0b0b0b", "#8a8984"
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.2, "font.size": 9,
                         "axes.titlesize": 10, "lines.linewidth": 2})
    powers = np.linspace(1, 200, 120)
    alone_curve = [R(P, "alone", pumps)["sq_right"] for P in powers]

    # ---------- Figure 1: ground rules ----------
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    phis = np.linspace(0, 2*np.pi, 241)
    reg = [solve_pump(0, 0, None, PUMP_MIRROR, ph)["used_in_ring"] for ph in phis]
    Dgrid = np.linspace(-1.5*PHC_BETA, 1.5*PHC_BETA, 301)
    phc = [max(solve_pump(PHC_BETA, D, None, PUMP_MIRROR, ph)["used_in_ring"] for D in Dgrid)
           for ph in phis]                    # laser re-tuned at each mirror phase
    ax[0].plot(phis/np.pi, reg, color=BLUE, label="regular ring")
    ax[0].plot(phis/np.pi, phc, color=ORANGE, lw=1.5, label="PhC ring (for contrast)")
    ax[0].set_xlabel("round-trip phase ring -> mirror -> ring  (units of $\\pi$)")
    ax[0].set_ylabel("fraction of input pump used in ring")
    ax[0].set_title("(a) Does the mirror distance matter?")
    ax[0].set_ylim(0, 1.1); ax[0].legend(frameon=False, loc="lower right")
    ax[0].text(0.05, 1.02, "regular: flat, CW and CCW just add", color=INK, fontsize=8)
    ax[0].text(0.05, 0.08, "PhC: ring reflects pump -> ring + mirror\nform a cavity, distance matters",
               color=INK, fontsize=8)

    rr = [R(P, "pumpR", pumps) for P in powers]
    oldb = [-10*np.log10(1 - eta*4*xx/(1+xx)**2)
            for xx in x_plain(powers)*(1 + 0.99*(1-2*eta)**2)]
    ax[1].plot(powers, alone_curve, color=MUTED, lw=1.5, label="ring alone (right)")
    ax[1].plot(powers, oldb, color=MUTED, lw=1.5, ls=":", label="old 'boost' model")
    ax[1].plot(powers, [r["sq_right"] for r in rr], color=BLUE, label="+ pump mirror: RIGHT")
    ax[1].plot(powers, [r["sq_left"] for r in rr], color=BLUE, ls="--", label="+ pump mirror: LEFT")
    ax[1].set_xlabel("input pump power (mW)"); ax[1].set_ylabel("squeezing out (dB)")
    ax[1].set_title("(b) Right and left squeezing are different")
    ax[1].legend(frameon=False, fontsize=8)

    losses = np.linspace(0, 0.2, 41)
    vals = []
    for L in losses:
        m = dict(PUMP_MIRROR); m["A_sq"], m["T_sq"] = L, 1-L
        vals.append(run(P0, "regular", right=m, pump=pumps["pumpR"])["sq_right"])
    ax[2].plot(losses*100, vals, color=BLUE)
    ax[2].axhline(R(P0, "alone", pumps)["sq_right"], color=MUTED, lw=1.5, ls=":", label="no mirror")
    ax[2].axvline(MIRROR_LOSS*100, color=MUTED, lw=1, ls="--")
    ax[2].text(MIRROR_LOSS*100+0.4, min(vals)+0.05, "placeholder\n(2%)", fontsize=8, color=INK)
    ax[2].set_xlabel("mirror loss in the squeezed bands (%)")
    ax[2].set_ylabel(f"squeezing out RIGHT at {P0:.0f} mW (dB)")
    ax[2].set_title("(c) Reflector loss costs right-side squeezing")
    ax[2].legend(frameon=False, fontsize=8)
    fig.suptitle("Ground rules: regular ring, pump from the left, pump mirror on the right",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_DIR/"ground_rules.png", dpi=140)

    # ---------- Figure 2: the three simulations ----------
    fig, ax = plt.subplots(2, 3, figsize=(14, 8.4))
    def curves(name):
        rows = [R(P, name, pumps, best_ph[name]) for P in powers]
        ok = np.array([max(r["x_ccw"], r["x_cw"]) < 1 for r in rows])
        return powers[ok], [r for r, k in zip(rows, ok) if k]

    # Sim 1
    for name, ring, col in (("alone", "regular", BLUE), ("sim1", "PhC", ORANGE)):
        Ps, rows = curves(name)
        ax[0, 0].plot(Ps, [r["sq_right"] for r in rows], color=col, label=f"{ring}: RIGHT")
        ax[0, 0].plot(Ps, [r["sq_left"] for r in rows], color=col, ls="--", label=f"{ring}: LEFT")
    ax[0, 0].set_title("Sim 1: PhC ring vs regular ring, no mirrors")
    # Sim 2 and 3 (regular ring)
    for j, name, title in ((1, "sim2", "Sim 2: pump mirror R + squeezing mirror L"),
                           (2, "sim3", "Sim 3: full mirror on the RIGHT")):
        Ps, rows = curves(name)
        ax[0, j].plot(powers, alone_curve, color=MUTED, lw=1.5, ls=":", label="ring alone (right)")
        ax[0, j].plot(Ps, [r["sq_right"] for r in rows], color=BLUE, label="RIGHT")
        ax[0, j].plot(Ps, [r["sq_left"] for r in rows], color=BLUE, ls="--", label="LEFT")
        ax[0, j].set_title(title + "  (best phase)")
    for j in range(3):
        ax[0, j].axhline(-10*np.log10(1-eta), color=MUTED, lw=1, ls=":")
        ax[0, j].set_xlabel("input pump power (mW)"); ax[0, j].set_ylabel("squeezing out (dB)")
        ax[0, j].set_ylim(-0.5, 6.5)
        ax[0, j].legend(frameon=False, fontsize=8, loc="upper left")
    ax[0, 0].text(195, -10*np.log10(1-eta)+0.12, "6 dB ceiling", ha="right", fontsize=8, color=INK)

    # (d) pump budget
    names = [("alone", "ring\nalone"), ("sim1", "Sim 1\nPhC"), ("sim2", "Sim 2"), ("sim3", "Sim 3")]
    xs = np.arange(len(names)); bot = np.zeros(len(names))
    for key, lab, colr, hatch in (("used_in_ring", "used in ring", "#3d3d3a", None),
                                  ("pump_left", "exits left", "#9a9994", None),
                                  ("pump_right", "exits right", "#d6d5d0", None),
                                  ("lost_in_mirrors", "lost in mirrors", "#ffffff", "///")):
        v = np.array([pumps[n][key] for n, _ in names])
        if hatch:
            ax[1, 0].bar(xs, v, bottom=bot, color=colr, edgecolor=MUTED, linewidth=0.8,
                         hatch=hatch, label=lab, width=0.6)
        else:
            ax[1, 0].bar(xs, v, bottom=bot, color=colr, edgecolor="white", linewidth=2,
                         label=lab, width=0.6)
        bot += v
    for k, (n, _) in enumerate(names):
        u = pumps[n]["used_in_ring"]
        ax[1, 0].text(k, u/2, f"{u:.0%}", ha="center", va="center", color="white", fontsize=8)
    ax[1, 0].set_xticks(xs); ax[1, 0].set_xticklabels([l for _, l in names], fontsize=8)
    ax[1, 0].set_ylabel("fraction of input pump"); ax[1, 0].set_ylim(0, 1.05)
    ax[1, 0].set_title("(d) Where the pump goes"); ax[1, 0].grid(False)
    ax[1, 0].legend(frameon=False, fontsize=7, ncol=2, loc="upper center",
                    bbox_to_anchor=(0.5, -0.14))

    # (e) phase sensitivity
    ph, v2 = phase_scan(P0, "sim2", pumps)
    ph, v3 = phase_scan(P0, "sim3", pumps)
    ax[1, 1].plot(ph/np.pi, v2, color=BLUE, label="Sim 2: RIGHT")
    ax[1, 1].plot(ph/np.pi, v3, color=BLUE, ls="--", label="Sim 3: LEFT")
    ax[1, 1].axhline(R(P0, "alone", pumps)["sq_right"], color=MUTED, lw=1.5, ls=":",
                     label="ring alone")
    ax[1, 1].set_xlabel("phase of the returning squeezed pair  (units of $\\pi$)")
    ax[1, 1].set_ylabel(f"squeezing at useful port, {P0:.0f} mW (dB)")
    ax[1, 1].set_title("(e) Phase of the returning squeezed pair")
    ax[1, 1].legend(frameon=False, fontsize=8)

    # (f) escape-efficiency what-if
    etas = np.linspace(0.55, 0.99, 23)
    s_alone, s2, s3 = [], [], []
    for e in etas:
        set_eta(e)
        pe = {n: best_pump(0.0, **LAYOUTS[n][1]) for n in ("alone", "sim2", "sim3")}
        s_alone.append(R(P0, "alone", pe)["sq_right"])
        s2.append(max(R(P0, "sim2", pe, q)["sq_right"] for q in np.linspace(0, 2*np.pi, 37)))
        s3.append(max(R(P0, "sim3", pe, q)["sq_left"] for q in np.linspace(0, 2*np.pi, 37)))
    set_eta(0.75)
    ax[1, 2].plot(etas, s_alone, color=MUTED, lw=1.5, ls=":", label="ring alone: RIGHT")
    ax[1, 2].plot(etas, s2, color=BLUE, label="Sim 2: RIGHT")
    ax[1, 2].plot(etas, s3, color=BLUE, ls="--", label="Sim 3: LEFT")
    ax[1, 2].axvline(0.75, color=MUTED, lw=1, ls="--")
    ax[1, 2].text(0.755, 0.3, "Vaidya ring", fontsize=8, color=INK)
    ax[1, 2].set_xlabel("escape efficiency $\\eta$ (same loaded Q)")
    ax[1, 2].set_ylabel(f"squeezing at useful port, {P0:.0f} mW (dB)")
    ax[1, 2].set_title("(f) Returning light must pass the ring: needs $\\eta \\to 1$")
    ax[1, 2].legend(frameon=False, fontsize=8)

    fig.suptitle(f"Three simulations, same Vaidya ring  (mirrors: R = {MIRROR_R:.0%}, "
                 f"loss = {MIRROR_LOSS:.0%};  PhC splitting = {2*PHC_BETA:.0f} linewidths)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_DIR/"three_simulations.png", dpi=140)
    print("Figures saved next to the script: ground_rules.png, three_simulations.png")

    # ---------- Figure 3: the five scenarios side by side ----------
    S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1.6, 1]})
    labels = [t.split(". ", 1)[1].replace(" + ", "\n+ ") for t, _ in summary]
    xs = np.arange(len(summary)); w = 0.2
    series = [("squeezing out LEFT",      "sq_left",  S1),
              ("made in ring, CCW",       "made_ccw", S2),
              ("made in ring, CW",        "made_cw",  S3),
              ("squeezing out RIGHT",     "sq_right", S4)]
    for k, (lab, key, col) in enumerate(series):
        v = np.array([max(r[key], 0) for _, r in summary])
        ax[0].bar(xs + (k-1.5)*w, v, w*0.9, color=col, label=lab)
        for xx, vv in zip(xs + (k-1.5)*w, v):
            if vv > 0.05:
                ax[0].text(xx, vv+0.1, f"{vv:.1f}", ha="center", fontsize=7, color=INK)
    ax[0].set_xticks(xs); ax[0].set_xticklabels(labels, fontsize=8)
    ax[0].set_ylabel("squeezing (dB below shot noise)")
    ax[0].set_title(f"Squeezing on each side and in the ring  ({P0:.0f} mW in)")
    ax[0].legend(frameon=False, fontsize=8, ncol=2); ax[0].set_ylim(0, 10)
    ax[0].grid(axis="x", visible=False)

    bot = np.zeros(len(summary))
    for key, lab, colr, hatch in (("used_in_ring", "used in ring", "#3d3d3a", None),
                                  ("pump_left", "exits left", "#9a9994", None),
                                  ("pump_right", "exits right", "#d6d5d0", None),
                                  ("lost_in_mirrors", "lost in mirrors", "#ffffff", "///")):
        v = np.array([P0*r["pump"][key] for _, r in summary])
        kw = dict(edgecolor=MUTED, linewidth=0.8, hatch=hatch) if hatch else \
             dict(edgecolor="white", linewidth=2)
        ax[1].bar(xs, v, bottom=bot, color=colr, label=lab, width=0.6, **kw)
        bot += v
    for k, (_, r) in enumerate(summary):
        u = P0*r["pump"]["used_in_ring"]
        ax[1].text(k, u/2, f"{u:.0f} mW", ha="center", va="center", color="white", fontsize=8)
    ax[1].set_xticks(xs); ax[1].set_xticklabels([f"{k+1}" for k in xs])
    ax[1].set_xlabel("scenario"); ax[1].set_ylabel("pump power (mW)")
    ax[1].set_title(f"Where the {P0:.0f} mW of pump goes"); ax[1].set_ylim(0, 105)
    ax[1].grid(False)
    ax[1].legend(frameon=False, fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    fig.tight_layout()
    fig.savefig(OUT_DIR/"five_scenarios.png", dpi=140)
    print("Comparison figure saved: five_scenarios.png")
