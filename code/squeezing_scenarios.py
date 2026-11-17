#!/usr/bin/env python3
# =====================================================================
#  squeezing_scenarios.py
#  Two-mode squeezing in a SiN microring: 8 layouts + the two-point bound
# =====================================================================
#
#  ONE FILE = Yaniv's engine (ring_three_configs.py, physics unchanged)
#           + the new PhC-ring-with-reflector layouts
#           + the fair-comparison tools (eta sweep, best eta, pump vs target).
#
#  The old "boost" model (ring_vs_reflector.py) is NOT included: it added
#  the recycled pump to the forward direction, which is not physical.
#
# ---------------------------------------------------------------------
#  HOW TO RUN
# ---------------------------------------------------------------------
#    python3 squeezing_scenarios.py              -> section A (seconds)
#    python3 squeezing_scenarios.py --eta        -> section B (seconds)
#    python3 squeezing_scenarios.py --power      -> section C (~1-2 min)
#    python3 squeezing_scenarios.py --fair       -> section D (~20-30 min)
#    python3 squeezing_scenarios.py --fair --quick   (coarse grid, ~3 min)
#    python3 squeezing_scenarios.py --all        -> A + B + C + D
#    add  --loss 0.005  to change the reflector loss (default 0.02)
#
# ---------------------------------------------------------------------
#  MAP OF THE FILE  (search for the "PART n" banners)
# ---------------------------------------------------------------------
#    PART 1  Fixed numbers: Vaidya ring, reflectors, PhC splitting
#    PART 2  Ring state that changes with eta  -> set_ring()
#    PART 3  Drive strength x                  -> x_plain()
#    PART 4  Classical pump: where the pump goes -> solve_pump(), best_pump()
#    PART 5  Quantum part: squeezer as a Bogoliubov map -> ring_squeezer()
#    PART 6  Reflectors                        -> make_mirror()
#    PART 7  THE LAYOUTS (all scenarios live here) -> LAYOUTS
#    PART 8  One run of one layout             -> run()
#    PART 9  Analysis tools: best phase, threshold, pump for a target,
#            best eta, two-point coupler bound
#    PART 10 Reports (the tables in the doc)   -> report_*()
#    PART 11 Command line                      -> main()
#
# ---------------------------------------------------------------------
#  THE LAYOUTS  (key in LAYOUTS -> what it is -> where the useful light exits)
# ---------------------------------------------------------------------
#    "plain"        Plain ring                               -> right
#    "ring_P"       Ring + pump reflector (right)            -> right
#    "phc"          PhC ring, no reflector                   -> right
#    "ring_P_S"     Ring + pump refl. (R) + squeeze refl. (L)-> right
#    "ring_F"       Ring + full reflector (right)            -> left
#    "phc_P"        PhC ring + pump reflector (right)        -> best of both
#    "phc_F"        PhC ring + full reflector (right)        -> left
#    "phc_P_S"      PhC ring + pump refl. (R) + squeeze (L)  -> right
#    two-point coupler: not a layout of the engine; an analytic bound,
#    see two_point_bound_dB() in PART 9.
#
#  Old names, for reading Yaniv's notes / the earlier chat:
#    alone=plain, pumpR=ring_P, sim1=phc, sim2=ring_P_S, sim3=ring_F,
#    aside=phc_F;  scenarios 6/7/8 (A/B/C) = phc_P / phc_F / phc_P_S
#
# ---------------------------------------------------------------------
#  CONVENTIONS
# ---------------------------------------------------------------------
#  * Pump enters from the LEFT, moving right on the bus.
#  * CCW mode = driven by light moving right; its pairs exit RIGHT.
#    CW  mode = driven by light moving left (from a right reflector);
#    its pairs exit LEFT.  ("direction rule")
#  * Squeezing is reported in dB below shot noise, at line centre
#    (Omega = 0), in the BUS waveguide (detection loss not included).
#  * Fair comparison: INTRINSIC Q fixed at 8e5 (Vaidya); the gap (eta)
#    is the designer's choice, the same eta for pump, signal and idler.
# =====================================================================

import argparse
import numpy as np
from scipy.optimize import minimize


# =====================================================================
# PART 1 - Fixed numbers
# =====================================================================
# Ring (Vaidya et al., Sci. Adv. 2020) - set by fabrication
Q_INTRINSIC = 8.0e5            # fixed in the fair comparison
ETA_VAIDYA  = 0.75             # escape efficiency of the published ring
wavelength  = 1550.2e-9
hbar, c     = 1.054e-34, 3.0e8
ring_radius, FSR = 120e-6, 190e9

w      = 2 * np.pi * c / wavelength
length = 2 * np.pi * ring_radius
n_g    = c / (FSR * length)
Lambda = hbar * w * (c / n_g) ** 2 * 1.0 / length   # nonlinear strength

# Reflectors - PLACEHOLDERS until there is a real grating design.
# In the reflect band:  R = 1 - 1.5 L,  T = 0.5 L,  A = L
#   (L = 2%  ->  R = 97%, T = 1%, A = 2%: Yaniv's placeholder)
# In the pass band:     T = 1 - L
MIRROR_LOSS_DEFAULT = 0.02

# PhC ring: splits the PUMP resonance only (signal/idler stay travelling).
# Splitting = 10 Vaidya linewidths (~10 GHz), FIXED IN ABSOLUTE TERMS:
# in units of the loaded linewidth it shrinks as eta (and the linewidth) grows.
PHC_SPLIT_VAIDYA_LINEWIDTHS = 10.0


# =====================================================================
# PART 2 - Ring state that changes with eta
# =====================================================================
# Everything below reads these globals. Call set_ring(eta) before any run.
eta = ETA_VAIDYA                  # escape efficiency (gap)
Q = Q_INTRINSIC * (1 - eta)       # loaded Q
half_linewidth = w / (2 * Q)      # loaded amplitude decay rate
PHC_BETA = PHC_SPLIT_VAIDYA_LINEWIDTHS / 2   # CW<->CCW coupling, in linewidths


def set_ring(new_eta, fixed="intrinsic"):
    """Set the escape efficiency (the coupling gap).
    fixed="intrinsic": intrinsic Q stays 8e5 (fair comparison, default).
    fixed="loaded":    loaded Q stays 2e5 (Yaniv's older 'what if')."""
    global eta, Q, half_linewidth, PHC_BETA
    eta = new_eta
    Q = Q_INTRINSIC * (1 - eta) if fixed == "intrinsic" else 2.0e5
    half_linewidth = w / (2 * Q)          # must follow Q (was fixed at import)
    # absolute splitting fixed -> in loaded linewidths it scales with Q
    PHC_BETA = (PHC_SPLIT_VAIDYA_LINEWIDTHS / 2) * Q / 2.0e5


# =====================================================================
# PART 3 - Drive strength x   (x = 1 is the threshold)
# =====================================================================
def x_plain(power_mW):
    """Drive x of a REGULAR ring with ALL pump in one direction, at the
    current eta and Q. Every layout scales this by its own stored pump:
        x_dir = x_plain * N_dir / N_ref   (PART 8)."""
    pump_inside = 4 * eta * Q * (power_mW / 1000) / (hbar * w ** 2)
    return Lambda * pump_inside / half_linewidth


# =====================================================================
# PART 4 - Classical pump: where does the pump go?
# =====================================================================
# Units: loaded linewidth kappa = 1, input amplitude s = 1 (|s|^2 = P_in).
# Unknowns: a1 = CCW pump, a2 = CW pump, SR = field hitting the ring from
# the left (moving right), SL = from the right (moving left).
#   (1/2 - iD) a1 - i beta a2 = sqrt(eta) SR      ring, CCW
#   (1/2 - iD) a2 - i beta a1 = sqrt(eta) SL      ring, CW
#   OR = SR - sqrt(eta) a1                        leaves to the right
#   OL = SL - sqrt(eta) a2                        leaves to the left
#   SR = tL * s                                   (left mirror passes pump)
#   SL = rR exp(i phiR) * OR                      (right mirror sends back)
# Regular ring = beta 0.  PhC ring = beta > 0 (standing waves, split).
def solve_pump(beta=0.0, detuning=0.0, left=None, right=None, phiR=0.0):
    TL = 1.0 if left is None else left["T_pump"]
    AL = 0.0 if left is None else left["A_pump"]
    RR = 0.0 if right is None else right["R_pump"]
    TR = 1.0 if right is None else right["T_pump"]
    AR = 0.0 if right is None else right["A_pump"]
    se, D, rRe = np.sqrt(eta), detuning, np.sqrt(RR) * np.exp(1j * phiR)
    tL = np.sqrt(TL)

    M = np.array([[0.5 - 1j * D, -1j * beta, -se, 0],
                  [-1j * beta, 0.5 - 1j * D, 0, -se],
                  [0, 0, 1, 0],
                  [rRe * se, 0, -rRe, 1]], dtype=complex)
    b = np.array([0, 0, tL, 0], dtype=complex)
    a1, a2, SR, SL = np.linalg.solve(M, b)
    OR, OL = SR - se * a1, SL - se * a2

    out = dict(a_ccw=a1, a_cw=a2,
               N_ccw=abs(a1) ** 2, N_cw=abs(a2) ** 2,       # stored pump
               pump_right=TR * abs(OR) ** 2,                 # exits right
               pump_left=TL * abs(OL) ** 2,                  # back to laser
               used_in_ring=(1 - eta) * (abs(a1) ** 2 + abs(a2) ** 2),
               lost_in_mirrors=AL * (1 + abs(OL) ** 2) + AR * abs(OR) ** 2)
    out["energy_check"] = (out["pump_right"] + out["pump_left"]
                           + out["used_in_ring"] + out["lost_in_mirrors"])
    return out


def best_pump(beta, left=None, right=None):
    """Choose the laser detuning D (and, for a PhC ring with a right
    reflector, the reflector phase phiR) that puts the MOST pump to use.
    Coarse grid, then Nelder-Mead polish.
    Regular ring: D = 0, and the reflector phase does not matter."""
    if beta == 0:
        s = solve_pump(0.0, 0.0, left, right, 0.0)
        s["D"], s["phiR"] = 0.0, 0.0
        return s
    scan = right is not None            # the phase is the coupling knob
    Ds = np.linspace(-1.5 * beta, 1.5 * beta, 241)
    phs = np.linspace(0, 2 * np.pi, 61) if scan else [0.0]
    best = None
    for D in Ds:
        for p in phs:
            u = solve_pump(beta, D, left, right, p)["used_in_ring"]
            if best is None or u > best[0]:
                best = (u, D, p)
    f = lambda v: -solve_pump(beta, v[0], left, right,
                              v[1] if scan else 0.0)["used_in_ring"]
    res = minimize(f, [best[1], best[2]], method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-10))
    D, p = res.x[0], (res.x[1] if scan else 0.0)
    s = solve_pump(beta, D, left, right, p)
    s["D"], s["phiR"] = D, p
    return s


# =====================================================================
# PART 5 - Quantum part: squeezed light as a Bogoliubov map
# =====================================================================
# Every field = sum_k alpha_k c_k + beta_k c_k^dagger, the c_k being
# independent vacuum inputs (bus, ring loss, mirror loss).
class Vacua:
    def __init__(self, n_max=60):
        self.n, self.n_max = 0, n_max

    def new(self):
        al, be = np.zeros(self.n_max, complex), np.zeros(self.n_max, complex)
        al[self.n] = 1.0
        self.n += 1
        return (al, be)


def dag(f):        return (np.conj(f[1]), np.conj(f[0]))
def scale(f, z):   return (z * f[0], z * f[1])
def add(*fs):      return (sum(f[0] for f in fs), sum(f[1] for f in fs))


def through_loss(f, keep, vac):
    """Keep a fraction 'keep' of the power; the rest becomes fresh vacuum
    (loss = beam splitter)."""
    return add(scale(f, np.sqrt(keep)), scale(vac.new(), np.sqrt(1 - keep)))


def ring_squeezer(fs, fi, x, theta, vac):
    """ONE direction of the ring at line centre (Collett-Gardiner).
    fs, fi : signal, idler ARRIVING on the bus (vacuum, or squeezed light
             sent back by a reflector - that is how cascades are built)
    x      : drive strength of this direction, theta : its pump phase
    Returns signal, idler LEAVING on the bus.
    For the quiet quadrature this is V_out = t^2 V_in + l^2 with
    t = 2 eta/(1+x) - 1  (t ~ 0 erases what came in)."""
    G = x / 2
    A = np.array([[0.5, -G * np.exp(1j * theta)],
                  [-G * np.exp(-1j * theta), 0.5]])
    Ainv = np.linalg.inv(A)
    ls, li = vac.new(), vac.new()                          # ring-loss vacua
    drive = [add(scale(fs, np.sqrt(eta)), scale(ls, np.sqrt(1 - eta))),
             add(scale(dag(fi), np.sqrt(eta)), scale(dag(li), np.sqrt(1 - eta)))]
    cav = [add(scale(drive[0], Ainv[r, 0]), scale(drive[1], Ainv[r, 1]))
           for r in range(2)]
    out_s = add(fs, scale(cav[0], -np.sqrt(eta)))
    out_id = add(dag(fi), scale(cav[1], -np.sqrt(eta)))
    return out_s, dag(out_id)


def squeezing_dB(bs, bi):
    """Best two-mode squeezing (dB below shot noise) and anti-squeezing:
    smallest / largest eigenvalue of the signal-idler covariance matrix
    = the best homodyne setting."""
    rows = []
    for b in (bs, bi):
        u, v = (b[0] + b[1]) / 2, 1j * (b[0] - b[1]) / 2
        rows.append(np.concatenate([2 * u.real, 2 * v.real]))   # q
        rows.append(np.concatenate([2 * u.imag, 2 * v.imag]))   # p
    R = np.array(rows)
    ev = np.linalg.eigvalsh(R @ R.T)
    return -10 * np.log10(ev[0]), 10 * np.log10(ev[-1])


# =====================================================================
# PART 6 - Reflectors
# =====================================================================
def make_mirror(reflects_pump, reflects_squeezing, loss=MIRROR_LOSS_DEFAULT):
    """Bragg reflector on the bus. Bands: 'pump' (1550 nm) and
    'sq' (1540/1560 nm). In its reflect band R = 1 - 1.5 L, T = 0.5 L;
    in a pass band R = 0, T = 1 - L. Loss A = L in both."""
    m = {}
    for band, refl in (("pump", reflects_pump), ("sq", reflects_squeezing)):
        R = (1 - 1.5 * loss) if refl else 0.0
        m["R_" + band], m["A_" + band] = R, loss
        m["T_" + band] = 1 - R - loss
    return m


# =====================================================================
# PART 7 - THE LAYOUTS   <<< every scenario is defined here >>>
# =====================================================================
# Each entry: (description, ring type, which reflectors, useful port)
#   reflectors: "P" = pump only, "S" = squeezing only, "F" = full
#   useful port "both" = two equivalent outputs, take the better one
LAYOUT_SPEC = {
    "plain":    ("Plain ring",                            "regular", {},                     "right"),
    "ring_P":   ("Ring + pump reflector",                 "regular", {"right": "P"},          "right"),
    "phc":      ("PhC ring, no reflector",                "PhC",     {},                     "right"),
    "ring_P_S": ("Ring + pump + squeeze reflectors",      "regular", {"left": "S", "right": "P"}, "right"),
    "ring_F":   ("Ring + full reflector",                 "regular", {"right": "F"},          "left"),
    "phc_P":    ("PhC ring + pump reflector",             "PhC",     {"right": "P"},          "both"),
    "phc_F":    ("PhC ring + full reflector",             "PhC",     {"right": "F"},          "left"),
    "phc_P_S":  ("PhC ring + pump & squeeze reflectors",  "PhC",     {"left": "S", "right": "P"}, "right"),
}
ORDER = list(LAYOUT_SPEC)


def layout(name, loss=MIRROR_LOSS_DEFAULT):
    """Turn a LAYOUT_SPEC entry into (ring type, mirror dicts, port)."""
    desc, ring, refl, port = LAYOUT_SPEC[name]
    kind = {"P": (True, False), "S": (False, True), "F": (True, True)}
    mirrors = {side: make_mirror(*kind[k], loss=loss) for side, k in refl.items()}
    return ring, mirrors, port


# =====================================================================
# PART 8 - One run of one layout (Yaniv's run(), physics unchanged)
# =====================================================================
def run(power_mW, ring="regular", left=None, right=None, sq_phase=0.0, pump=None):
    """Pump bookkeeping + squeezing on each side for one layout.
    sq_phase: extra round-trip phase of a RETURNING squeezed pair relative
              to the pump (0 = 'in step'). Only matters when a reflector
              sends squeezed light back (ring_P_S, ring_F, phc_F, phc_P_S).
    Three branches:
      (1) squeeze reflector on the LEFT  -> CW pairs go left, bounce, and
          are fed into the CCW squeezer; output RIGHT   (ring_P_S, phc_P_S)
      (2) full reflector on the RIGHT    -> CCW pairs go right, bounce, and
          are fed into the CW squeezer; output LEFT     (ring_F, phc_F)
      (3) otherwise: two independent squeezers, CCW -> right, CW -> left
                                                (plain, ring_P, phc, phc_P)"""
    beta = PHC_BETA if ring == "PhC" else 0.0
    if pump is None:
        pump = best_pump(beta, left, right)

    # drive of each direction, scaled from the plain ring by stored pump
    x0 = x_plain(power_mW)
    N_ref = 4 * eta                       # plain ring, all pump one way
    x_ccw = x0 * pump["N_ccw"] / N_ref
    x_cw = x0 * pump["N_cw"] / N_ref
    th_ccw, th_cw = 2 * np.angle(pump["a_ccw"]), 2 * np.angle(pump["a_cw"])

    vac = Vacua()
    fs, fi = vac.new(), vac.new()         # vacuum arriving from the left
    if left is not None:                  # ... through the left mirror
        fs, fi = through_loss(fs, left["T_sq"], vac), through_loss(fi, left["T_sq"], vac)

    if left is not None and left["R_sq"] > 0:
        # ---- branch (1): squeeze reflector on the left -------------------
        gs, gi = vac.new(), vac.new()                     # vacuum from right
        cw_s, cw_i = ring_squeezer(gs, gi, x_cw, th_cw, vac)   # weak CW source
        back_s = through_loss(cw_s, left["R_sq"], vac)          # bounce
        back_i = through_loss(cw_i, left["R_sq"], vac)
        ph = np.exp(1j * (sq_phase / 2))                  # free phase (heater)
        back_s, back_i = scale(back_s, ph), scale(back_i, ph)
        fs = add(back_s, scale(vac.new(), np.sqrt(left["T_sq"])))
        fi = add(back_i, scale(vac.new(), np.sqrt(left["T_sq"])))
        ccw_s, ccw_i = ring_squeezer(fs, fi, x_ccw, th_ccw, vac)  # crosses CCW
        right_s, right_i = ccw_s, ccw_i
        if right is not None:                             # pump mirror loss
            right_s = through_loss(ccw_s, right["T_sq"], vac)
            right_i = through_loss(ccw_i, right["T_sq"], vac)
        left_s = through_loss(cw_s, left["T_sq"], vac)
        left_i = through_loss(cw_i, left["T_sq"], vac)
    else:
        ccw_s, ccw_i = ring_squeezer(fs, fi, x_ccw, th_ccw, vac)  # CCW source
        if right is not None and right["R_sq"] > 0:
            # ---- branch (2): full reflector on the right -----------------
            # pump took round-trip phase phiR; the pair takes 2 phiR - sq_phase
            ph = np.exp(1j * (pump["phiR"] - sq_phase / 2))
            back_s = scale(through_loss(ccw_s, right["R_sq"], vac), ph)
            back_i = scale(through_loss(ccw_i, right["R_sq"], vac), ph)
            right_s = through_loss(ccw_s, right["T_sq"], vac)
            right_i = through_loss(ccw_i, right["T_sq"], vac)
        else:
            # ---- branch (3): nothing comes back --------------------------
            back_s, back_i = vac.new(), vac.new()
            if right is not None:          # pump-only mirror: pass with loss
                right_s = through_loss(ccw_s, right["T_sq"], vac)
                right_i = through_loss(ccw_i, right["T_sq"], vac)
            else:
                right_s, right_i = ccw_s, ccw_i
        cw_s, cw_i = ring_squeezer(back_s, back_i, x_cw, th_cw, vac)  # CW
        left_s, left_i = cw_s, cw_i
        if left is not None:
            left_s = through_loss(left_s, left["T_sq"], vac)
            left_i = through_loss(left_i, left["T_sq"], vac)

    sqR, antiR = squeezing_dB(right_s, right_i)
    sqL, antiL = squeezing_dB(left_s, left_i)
    return dict(pump=pump, x_ccw=x_ccw, x_cw=x_cw,
                sq_right=sqR, anti_right=antiR, sq_left=sqL, anti_left=antiL)


# =====================================================================
# PART 9 - Analysis tools
# =====================================================================
_pump_cache = {}


def pump_for(name, loss):
    """Pump solution for a layout at the CURRENT eta (cached: the pump is
    linear, so it does not depend on power)."""
    key = (name, round(eta, 6), loss, round(PHC_BETA, 9))
    if key not in _pump_cache:
        ring, mir, _ = layout(name, loss)
        beta = PHC_BETA if ring == "PhC" else 0.0
        _pump_cache[key] = best_pump(beta, **mir)
    return _pump_cache[key]


def returns_squeezing(name):
    return any(k in ("S", "F") for k in LAYOUT_SPEC[name][2].values())


def port_dB(r, port):
    if port == "both":
        return max(r["sq_left"], r["sq_right"])
    return r["sq_" + port]


PHASES = np.linspace(0, 2 * np.pi, 25)[:-1]


def squeeze_at(power_mW, name, loss=MIRROR_LOSS_DEFAULT, mode="best"):
    """Squeezing (dB) at the useful port, at the current eta.
    mode = "best" | "worst" | "instep" (phase of a returning pair).
    Returns nan above threshold (x >= 1 in either direction)."""
    ring, mir, port = layout(name, loss)
    pump = pump_for(name, loss)
    phases = PHASES if (returns_squeezing(name) and mode != "instep") else [0.0]
    vals = []
    for ph in phases:
        r = run(power_mW, ring, sq_phase=ph, pump=pump, **mir)
        vals.append(np.nan if max(r["x_ccw"], r["x_cw"]) >= 1 else port_dB(r, port))
    vals = np.array(vals)
    if np.all(np.isnan(vals)):
        return np.nan
    return np.nanmin(vals) if mode == "worst" else np.nanmax(vals)


def threshold_mW(name, loss=MIRROR_LOSS_DEFAULT):
    """Pump at which the stronger direction reaches x = 1 (current eta)."""
    ring, mir, _ = layout(name, loss)
    r = run(100.0, ring, pump=pump_for(name, loss), **mir)
    return 100.0 / max(r["x_ccw"], r["x_cw"])


def power_for_target(target_dB, name, loss=MIRROR_LOSS_DEFAULT, mode="best"):
    """Least pump (mW) that reaches target_dB at the current eta, below
    threshold. nan if the target cannot be reached at this eta.
    Bisection on log(P): squeezing rises monotonically up to threshold."""
    Pth = threshold_mW(name, loss)
    if not squeeze_at(Pth * 0.9999, name, loss, mode) >= target_dB:
        return np.nan
    lo, hi = 1e-6, Pth * 0.9999
    for _ in range(45):
        mid = np.sqrt(lo * hi)
        if squeeze_at(mid, name, loss, mode) >= target_dB:
            hi = mid
        else:
            lo = mid
    return hi


def best_eta_for_target(target_dB, name, loss=MIRROR_LOSS_DEFAULT, mode="best",
                        coarse=0.01, fine=0.001):
    """FAIR COMPARISON: scan eta (coarse, then fine around the best) and
    return (least pump in mW, the eta that gives it)."""
    def scan(etas):
        best = (np.inf, None)
        for e in etas:
            set_ring(e)
            P = power_for_target(target_dB, name, loss, mode)
            if np.isfinite(P) and P < best[0]:
                best = (P, e)
        return best
    P, e = scan(np.round(np.arange(0.50, 0.981, coarse), 4))
    if e is None:
        return np.nan, None
    if fine:
        P2, e2 = scan(np.round(np.arange(max(0.5, e - 1.2 * coarse),
                                         min(0.99, e + 1.2 * coarse) + 1e-9, fine), 4))
        if P2 < P:
            P, e = P2, e2
    set_ring(ETA_VAIDYA)
    return P, e


def best_eta_for_power(power_mW, name, loss=MIRROR_LOSS_DEFAULT, step=0.01):
    """Highest squeezing (dB) at a fixed pump power, over eta."""
    best = (-np.inf, None)
    for e in np.round(np.arange(0.50, 0.981, step), 4):
        set_ring(e)
        v = squeeze_at(power_mW, name, loss)
        if np.isfinite(v) and v > best[0]:
            best = (v, e)
    set_ring(ETA_VAIDYA)
    return best


# ---- Two-point (interferometric / MZI-coupled) coupler: ideal bound ----
# All pump absorbed (critical coupling for the pump), all of it in ONE
# direction, signal/idler escape with the gap's eta. No bypass loss.
#   x_bound = x_plain / (4 eta (1 - eta))
def two_point_bound_dB(power_mW, e):
    set_ring(e)
    x = x_plain(power_mW) / (4 * e * (1 - e))
    set_ring(ETA_VAIDYA)
    if x >= 1:
        return np.nan
    return -10 * np.log10(1 - 4 * e * x / (1 + x) ** 2)


def two_point_power_for_target(target_dB, e):
    """Closed form: solve 1 - 4 e x/(1+x)^2 = V for x < 1, then P."""
    V = 10 ** (-target_dB / 10)
    a, b, cc = (1 - V), 2 * (1 - V) - 4 * e, (1 - V)
    disc = b * b - 4 * a * cc
    if disc < 0:
        return np.nan
    x = (-b - np.sqrt(disc)) / (2 * a)
    if not 0 < x < 1:
        return np.nan
    set_ring(e)
    x100 = x_plain(100.0) / (4 * e * (1 - e))
    set_ring(ETA_VAIDYA)
    return 100.0 * x / x100


def two_point_best_eta_for_target(target_dB):
    best = (np.inf, None)
    for e in np.arange(0.5, 0.99, 0.001):
        P = two_point_power_for_target(target_dB, e)
        if np.isfinite(P) and P < best[0]:
            best = (P, round(e, 3))
    return best


# =====================================================================
# PART 10 - Reports (the tables of the comparison doc)
# =====================================================================
def report_operating_point(loss):
    """A. Vaidya's ring (eta 0.75), 100 mW: squeezing per port + pump budget."""
    set_ring(ETA_VAIDYA)
    print("=" * 110)
    print(f"A. At Vaidya's ring: eta = 0.75, 100 mW, reflector loss {100*loss:.1f}%")
    print("=" * 110)
    print(f"{'layout':40s}| {'L (dB)':>7} {'R (dB)':>7} | {'used':>6} {'out L':>6} {'out R':>6} {'lost':>5} (mW) | x CCW / CW")
    for name in ORDER:
        ring, mir, port = layout(name, loss)
        pump = pump_for(name, loss)
        if returns_squeezing(name):        # best phase of the returning pair
            vals = [(port_dB(run(100, ring, sq_phase=ph, pump=pump, **mir), port), ph) for ph in PHASES]
            ph = max(vals)[1]
        else:
            ph = 0.0
        r = run(100.0, ring, sq_phase=ph, pump=pump, **mir)
        p = r["pump"]
        print(f"{LAYOUT_SPEC[name][0]:40s}| {r['sq_left']:7.2f} {r['sq_right']:7.2f} | "
              f"{100*p['used_in_ring']:6.1f} {100*p['pump_left']:6.1f} {100*p['pump_right']:6.1f} "
              f"{100*p['lost_in_mirrors']:5.1f}      | {r['x_ccw']:.2f} / {r['x_cw']:.2f}")
    print(f"{'Two-point coupler (ideal bound)':40s}|    0.00 {two_point_bound_dB(100, 0.75):7.2f} | "
          f"  100.0    0.0    0.0   0.0      | {x_plain(100)/0.75:.2f} / 0.00")
    print("(useful port: see the LAYOUTS map at the top; returning pairs at their best phase)\n")


def report_eta_sweep(loss, etas=(0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)):
    """B. Same power (100 mW), different eta: dB at the useful port."""
    print("=" * 110)
    print(f"B. 100 mW, different eta, reflector loss {100*loss:.1f}%  (dB at the useful port)")
    print("=" * 110)
    print(f"{'eta':>5} | " + " | ".join(f"{n:>8}" for n in ORDER) + " | two-point")
    for e in etas:
        set_ring(e)
        row = []
        for n in ORDER:
            v = squeeze_at(100.0, n, loss)
            row.append("     OPO" if not np.isfinite(v) else f"{v:8.2f}")
        b = two_point_bound_dB(100.0, e)
        print(f"{e:5.2f} | " + " | ".join(row) + f" | {b:8.2f}")
    set_ring(ETA_VAIDYA)
    print("(OPO = above threshold at this eta)\n")


def report_best_eta_vs_power(loss, powers=(50, 100, 200, 400, 800)):
    """C. Same power, each layout at its best eta: dB (best eta)."""
    print("=" * 110)
    print(f"C. Best eta for each pump power, reflector loss {100*loss:.1f}%   cells: dB (eta)")
    print("=" * 110)
    print(f"{'P (mW)':>6} | " + " | ".join(f"{n:>12}" for n in ORDER) + " | two-point")
    for P in powers:
        row = []
        for n in ORDER:
            v, e = best_eta_for_power(P, n, loss)
            row.append(f"{v:5.2f} ({e:.2f})" if e is not None else "           -")
        bb = max(((two_point_bound_dB(P, e), e) for e in np.arange(0.5, 0.99, 0.002)),
                 key=lambda t: -np.inf if not np.isfinite(t[0]) else t[0])
        print(f"{P:6d} | " + " | ".join(row) + f" | {bb[0]:5.2f} ({bb[1]:.2f})")
    print()


def report_fair(loss, targets=(3, 4, 6, 8, 10, 12), quick=False):
    """D. FAIR COMPARISON: pump needed for a squeezing target, each layout
    at its best eta, relative to the plain ring at its best eta."""
    fine = 0 if quick else 0.001
    coarse = 0.02 if quick else 0.01
    print("=" * 110)
    print(f"D. Fair comparison (intrinsic Q fixed), reflector loss {100*loss:.1f}%")
    print("   cells: pump relative to the plain ring (best eta of that layout)")
    print("=" * 110)
    ref = {t: best_eta_for_target(t, "plain", loss, coarse=coarse, fine=fine) for t in targets}
    print("plain ring: " + "  ".join(f"{t} dB: {ref[t][0]:.0f} mW (eta {ref[t][1]})" for t in targets))
    for n in ORDER[1:]:
        modes = ["best", "instep"] if n == "phc_F" else ["best"]
        for mode in modes:
            cells = []
            for t in targets:
                P, e = best_eta_for_target(t, n, loss, mode, coarse=coarse, fine=fine)
                cells.append(f"{t}dB {P/ref[t][0]:.2f}x" if np.isfinite(P) else f"{t}dB  -  ")
            label = LAYOUT_SPEC[n][0] + (" [in-step]" if mode == "instep" else "")
            print(f"{label:46s}: " + " | ".join(cells), flush=True)
    cells = []
    for t in targets:
        P, e = two_point_best_eta_for_target(t)
        cells.append(f"{t}dB {P/ref[t][0]:.2f}x")
    print(f"{'Two-point coupler (ideal bound)':46s}: " + " | ".join(cells))
    print()


# =====================================================================
# PART 11 - Command line
# =====================================================================
def main():
    ap = argparse.ArgumentParser(description="Squeezing layouts: comparison tables")
    ap.add_argument("--eta", action="store_true", help="B: eta sweep at 100 mW")
    ap.add_argument("--power", action="store_true", help="C: best eta vs pump power")
    ap.add_argument("--fair", action="store_true", help="D: pump vs squeezing target")
    ap.add_argument("--all", action="store_true", help="A + B + C + D")
    ap.add_argument("--quick", action="store_true", help="coarse eta grid for D")
    ap.add_argument("--loss", type=float, default=MIRROR_LOSS_DEFAULT, help="reflector loss")
    a = ap.parse_args()

    report_operating_point(a.loss)
    if a.eta or a.all:
        report_eta_sweep(a.loss)
    if a.power or a.all:
        report_best_eta_vs_power(a.loss)
    if a.fair or a.all:
        report_fair(a.loss, quick=a.quick)


if __name__ == "__main__":
    main()
