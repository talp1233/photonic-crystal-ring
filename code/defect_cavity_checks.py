"""
defect_cavity_checks.py  --  "smaller mode volume" checks for the SiN squeezing ring
====================================================================================

Question checked
----------------
In the engine the drive parameter scales as  x ~ U / (V_eff * kappa).
Can a photonic-crystal (PhC) DEFECT cavity inside the ring shrink the mode volume
V_eff by a factor F (relative to the plain ring of Vaidya et al.) and so cut the
pump power needed for a squeezing target?

Three checks (each is one subcommand):

  gap     MPB, full 3D unit cell of the corrugated Vaidya waveguide
          (1.0 x 0.8 um SiN in SiO2).  Band gap vs corrugation depth, convergence
          between two resolutions, and the period retune that puts mid-gap at 1550 nm.

  cmt     Coupled-mode theory (numpy only, no MEEP needed).  Defect length that
          traps pump + signal + idler, decay lengths, and the nonlinear overlap of
          the three STANDING-wave modes -> realistic F.

  cavity  MEEP, 2D effective-index model of the defect cavity (grating - defect -
          grating, abrupt or tapered).  Mode frequencies and radiation Q (Harminv),
          then the frequency mismatch  f_s + f_i - 2 f_p.

  field   MEEP, same tapered cavity driven narrowband at ONE mode frequency; the
          DFT after the source is the mode profile (saved as mode_<tag>.npz).

  overlap F from the simulated fields: FWM overlap of idler/pump/signal profiles
          relative to a travelling ring mode with the same transverse profile.

  summary Reads the JSON outputs of the three checks and prints the tables and the
          net estimate  gain in x ~ F / (kappa_total / kappa_0).  No MEEP needed.

Results obtained so far (Nov 15, 2026)
--------------------------------------
  gap:     converged to < 1 % (res 32 vs 48);  ~1.1 THz gap at 100 nm, ~3.7 THz at
           200 nm corrugation; period must grow ~0.2 % (~1 nm) at 200 nm.
  cmt:     L_defect ~ 44 um, decay 12 um (pump) / 15 um (signal); length-only
           estimate F ~ 13.  The CMT overlap is very sensitive to sub-um details
           of the defect length (F swings 0.2 - 4.2 over 0.5 um), so it is only a guide.
  field/overlap (MEEP fields, tapered cavity): F = 3.2 - 3.4 (Ey only),
           2.4 - 2.5 (full vector) -> F ~ 3.
  cavity:  abrupt   Q_rad 3e4 - 7e4;  tapered (20 periods) Q_rad 1e5 - 2e6, and it
           differs between pump, signal and idler; mismatch 36 - 106 GHz
           vs ~1 GHz linewidth (defect not tuned).
  net:     required pump  P / P_ring ~ r_p * r_s / F,  r = 1 + Q_int / Q_rad of
           that mode.  With today's tapered Q_rad (some modes ~1e5 - 3e5) this
           is ~0.5 - 4, i.e. no gain.  Gain ~2x needs Q_rad >~ 4e6 on ALL three
           modes (r <~ 1.2) AND f_s + f_i - 2 f_p within ~1 GHz.
           Caveats: 2D (no vertical radiation), untuned defect.

How to run
----------
  MPB/MEEP live in the conda-forge package "pymeep" (Linux / macOS / WSL2):
      conda create -n mp -c conda-forge python=3.11 pymeep numpy
      conda activate mp

  python defect_cavity_checks.py gap                       # ~1 h on 4 cores
  python defect_cavity_checks.py cmt                       # seconds
  python defect_cavity_checks.py cavity --ntap 0  --tag abrupt   # ~20-40 min each
  python defect_cavity_checks.py cavity --ntap 20 --tag taper
  python defect_cavity_checks.py field 0.638609 m1            # ~20 min each, one per mode
  python defect_cavity_checks.py field 0.642288 m2            #   (frequencies from cav_taper.json)
  python defect_cavity_checks.py field 0.646087 m3
  python defect_cavity_checks.py field 0.649533 m4
  python defect_cavity_checks.py overlap m1 m2 m3 m4          # numpy only
  python defect_cavity_checks.py summary                   # tables from the JSON files

Units: lengths in um, frequencies in c/um (MEEP/MPB units).  f = 1/lambda[um].
1 c/um = 299.79 THz.
"""
import argparse, json, os, sys
import numpy as np

# ----------------------------------------------------------------------------
# PART 0. Common constants  (Vaidya-like ring, 1550 nm)
# ----------------------------------------------------------------------------
C_THZ   = 299.792458        # 1 c/um in THz
LAM0    = 1.55              # pump wavelength (um)
F0      = 1 / LAM0          # pump frequency (c/um)
N_SIN   = 1.996             # SiN index at 1550
N_OX    = 1.444             # SiO2 index
W0, T0  = 1.0, 0.8          # waveguide width and thickness (um)
N_EFF   = 1.774             # MPB cross-section result (xsec.py)
N_G     = 2.104             # MPB group index
N_2D    = 1.8752            # 2D effective slab index that reproduces N_EFF
R_RING  = 120.0             # ring radius (um)
L_RING  = 2 * np.pi * R_RING
Q_INT   = 8e5               # intrinsic Q used in all engine comparisons
DF_SIG  = 1.25              # signal/idler offset from the pump (THz), 1540/1560 nm
LINEWIDTH_GHZ = C_THZ * 1e3 * F0 / 2e5   # loaded linewidth of Vaidya's ring (~0.97 GHz)


def save(obj, name):
    json.dump(obj, open(name, "w"), indent=1)
    print("saved", name)


# ----------------------------------------------------------------------------
# PART 1. Band gap of the corrugated waveguide (MPB, 3D)
# ----------------------------------------------------------------------------
# Width modulation  w(x) = W0 + 2*delta*cos(2 pi x / a)   (both walls move by delta).
# The unit cell is built from ns thin blocks.  At k = 0.5 (Bragg point) the two
# lowest z-even bands are the band edges; their difference is the gap.
def mpb_gap(delta, a, res_yz, res_x, ns=64, nb=4):
    import meep as mp
    from meep import mpb
    geo = []
    for i in range(ns):
        x = (i + 0.5) / ns * a - a / 2
        w = W0 + 2 * delta * np.cos(2 * np.pi * x / a)
        geo.append(mp.Block(center=mp.Vector3(x, 0, 0), size=mp.Vector3(a / ns, w, T0),
                            material=mp.Medium(index=N_SIN)))
    ms = mpb.ModeSolver(geometry_lattice=mp.Lattice(size=mp.Vector3(a, 3.0, 3.0)),
                        geometry=geo, default_material=mp.Medium(index=N_OX),
                        resolution=mp.Vector3(res_x, res_yz, res_yz),
                        num_bands=nb, k_points=[mp.Vector3(0.5, 0, 0)])
    ms.run_zeven()                       # TE-like (Ex-dominant) modes
    return [float(v) for v in ms.get_freqs()]


def check_gap(args):
    out = {}
    a = LAM0 / (2 * N_EFF)               # first-order Bragg period, 0.43696 um
    for delta in args.deltas:
        for res in ((32, 160), (48, 240)):           # (transverse, along x) px/um
            f = mpb_gap(delta, a, *res)
            out[f"{delta}|{a:.5f}|{res}"] = f
            print(f"delta={delta} res={res}: gap = {(f[1]-f[0])*C_THZ:.3f} THz, "
                  f"mid = {0.5*(f[0]+f[1]):.5f} (target {F0:.5f})", flush=True)
        # retune: frequencies scale as 1/a, so a_new = a * f_mid / f0
        a2 = a * 0.5 * (f[0] + f[1]) / F0
        f = mpb_gap(delta, a2, 48, 240)
        out[f"{delta}|{a2:.5f}|tuned"] = f
        print(f"delta={delta} tuned a={a2:.5f}: mid = {0.5*(f[0]+f[1]):.5f}", flush=True)
    save(out, "check1.json")


# ----------------------------------------------------------------------------
# PART 2. Coupled-mode estimate of the defect cavity and of F
# ----------------------------------------------------------------------------
# kappa_b = pi n_g dnu_gap / c  (grating coupling, 1/um).  A mode detuned by delta
# from mid-gap decays in the grating with  gamma = sqrt(kappa_b^2 - delta^2).
# Choose the defect length so that the pump sits at mid-gap and the signal /
# idler at +-1.25 THz are the next resonances (round-trip phase = 2 pi).
# F = |overlap of the three normalized modes| * L_ring : the FWM coupling of the
# defect relative to a travelling-wave ring mode (overlap = 1/L_ring).
def check_cmt(args):
    dfs = DF_SIG * 1e12 / 3e14                   # signal offset in c/um (c = 3e8 m/s, as in the saved run;
                                                 # the result is sensitive to this at the 0.1 % level, see notes)
    out = {}
    for label, gap in [("2D", 0.0129), ("3D res32", 0.01244)]:     # gaps at 200 nm
        kap = np.pi * N_G * gap
        ds = 2 * np.pi * N_G * dfs
        u = ds / kap                             # signal position inside the half gap
        if u >= 1:
            out[label] = "signal outside the gap"; continue
        Ld = (np.pi - np.arcsin(u)) / ds
        gp, gs = kap, np.sqrt(kap**2 - ds**2)
        beta_B = 2 * np.pi * N_EFF * F0
        z = np.linspace(-Ld/2 - 12/gs, Ld/2 + 12/gs, 400001); dz = z[1] - z[0]

        def mode(delta, gam, parity):            # standing wave inside, decaying outside
            b = beta_B + delta
            inside = np.abs(z) <= Ld / 2
            car_in = np.cos(b * z) if parity == 0 else np.sin(b * z)
            zi = Ld / 2
            val = np.cos(b * zi) if parity == 0 else np.sin(b * zi)
            env = np.exp(-gam * (np.abs(z) - zi))
            sgn = np.where(z > 0, 1, 1 if parity == 0 else -1)
            car_out = sgn * np.cos(beta_B * (np.abs(z) - zi) + np.arccos(np.clip(val, -1, 1)))
            E = np.where(inside, car_in, env * car_out)
            return E / np.sqrt(np.sum(E**2) * dz)

        Ep, Es, Ei = mode(0.0, gp, 0), mode(+ds, gs, 1), mode(-ds, gs, 1)
        O = np.sum(Ep**2 * Es * Ei) * dz
        out[label] = dict(kappa_per_um=kap, decay_len_pump_um=1/gp, decay_len_signal_um=1/gs,
                          signal_over_halfgap=u, defect_len_um=Ld,
                          Leff_pump_um=Ld + 1/gp, Leff_signal_um=Ld + 1/gs,
                          F_vs_ring=abs(O) * L_RING, F_length_only=L_RING / (Ld + 1/gp))
    print(json.dumps(out, indent=1))
    save(out, "cmt.json")


# ----------------------------------------------------------------------------
# PART 3. Defect cavity in MEEP (2D effective index) -> frequencies and Q_rad
# ----------------------------------------------------------------------------
# Geometry along x:  PML | grating (Lg) | defect (Ld, plain width W0) | grating | PML
# ntap > 0 ramps the corrugation linearly over ntap periods next to the defect.
# Ey source with odd mirror symmetry in y (TE fundamental).  Harminv at two points
# returns the resonances and their Q (radiation + leakage through the grating).
def check_cavity(args):
    import meep as mp
    delta, Ld, Lg, ntap, res, T, a = args.delta, args.Ld, args.Lg, args.ntap, args.res, args.T, args.a
    nper = int(round(Lg / a))
    dpml, sy = 1.5, 5.0
    sx = Ld + 2 * nper * a + 2 * dpml + 3.0

    def width(x):
        s = abs(x) - Ld / 2                      # distance into the grating
        if s <= 0 or s >= nper * a:
            return W0
        k = int(s // a); u = (s - k * a) / a
        A = delta * min(1.0, (k + 0.5) / ntap) if ntap > 0 else delta
        return W0 + 2 * A * np.cos(2 * np.pi * u)

    gr, gy = 80, 2.4                             # material grid: 80 px/um, 2.4 um tall
    Nx, Ny = int(sx * gr), int(gy * gr)
    xs = (np.arange(Nx) + 0.5) / gr - sx / 2
    ys = (np.arange(Ny) + 0.5) / gr - gy / 2
    wv = np.array([width(x) for x in xs])
    weights = (np.abs(ys)[None, :] < wv[:, None] / 2).astype(float)
    mg = mp.MaterialGrid(mp.Vector3(Nx, Ny), mp.Medium(index=N_OX), mp.Medium(index=N_2D),
                         weights=weights, do_averaging=True, beta=0)
    geo = [mp.Block(center=mp.Vector3(), size=mp.Vector3(sx, gy), material=mg)]
    src = [mp.Source(mp.GaussianSource(F0, fwidth=0.016), component=mp.Ey, center=mp.Vector3(0.37, 0))]
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geo, sources=src, resolution=res,
                        boundary_layers=[mp.PML(dpml)], default_material=mp.Medium(index=N_OX),
                        symmetries=[mp.Mirror(mp.Y, phase=-1)])
    hs = [mp.Harminv(mp.Ey, p, F0, 0.016) for p in (mp.Vector3(1.1, 0), mp.Vector3(-7.3, 0))]
    sim.run(*[mp.after_sources(h) for h in hs], until_after_sources=T)
    modes = [[dict(f=m.freq, Q=m.Q, amp=abs(m.amp)) for m in h.modes] for h in hs]
    save(dict(delta=delta, Ld=Ld, Lg=Lg, ntap=ntap, res=res, T=T, a=a, modes=modes),
         f"cav_{args.tag}.json")


# ----------------------------------------------------------------------------
# PART 3b. Mode profiles (MEEP) and F from the simulated fields
# ----------------------------------------------------------------------------
# Same tapered cavity as PART 3.  A narrow Gaussian source (fwidth 0.0015) at the
# mode frequency excites mainly that mode; a DFT accumulated AFTER the source
# filters the rest.  Two sources (x = 0.37 and -5.13 um) so that both even and
# odd modes are excited.
def check_field(args):
    import meep as mp
    fm, tag = args.fm, args.tag
    delta, Ld, Lg, ntap, res, a = 0.2, 44.4, 80.0, 20, 25, 0.43987
    nper = int(round(Lg / a)); dpml, sy = 1.5, 5.0
    sx = Ld + 2 * nper * a + 2 * dpml + 3.0
    def width(x):
        s = abs(x) - Ld / 2
        if s <= 0 or s >= nper * a: return W0
        k = int(s // a); u = (s - k * a) / a
        return W0 + 2 * delta * min(1.0, (k + 0.5) / ntap) * np.cos(2 * np.pi * u)
    gr, gy = 80, 2.4; Nx, Ny = int(sx * gr), int(gy * gr)
    xs = (np.arange(Nx) + 0.5) / gr - sx / 2; ys = (np.arange(Ny) + 0.5) / gr - gy / 2
    wv = np.array([width(x) for x in xs]); wts = (np.abs(ys)[None, :] < wv[:, None] / 2).astype(float)
    mg = mp.MaterialGrid(mp.Vector3(Nx, Ny), mp.Medium(index=N_OX), mp.Medium(index=N_2D),
                         weights=wts, do_averaging=True, beta=0)
    geo = [mp.Block(center=mp.Vector3(), size=mp.Vector3(sx, gy), material=mg)]
    fw = 0.0015
    src = [mp.Source(mp.GaussianSource(fm, fwidth=fw), component=mp.Ey, center=mp.Vector3(x0, 0))
           for x0 in (0.37, -5.13)]
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geo, sources=src, resolution=res,
                        boundary_layers=[mp.PML(dpml)], default_material=mp.Medium(index=N_OX),
                        symmetries=[mp.Mirror(mp.Y, phase=-1)])
    vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(Ld + 2 * nper * a, gy))
    sim.run(until=src[0].src.swigobj.last_time_max())          # until the source is off
    dft = sim.add_dft_fields([mp.Ex, mp.Ey], fm, 0, 1, where=vol)
    sim.run(until=4000)                                          # ring-down, DFT at fm
    np.savez(f"mode_{tag}.npz", Ex=sim.get_dft_array(dft, mp.Ex, 0), Ey=sim.get_dft_array(dft, mp.Ey, 0),
             eps=sim.get_array(vol=vol, component=mp.Dielectric), fm=fm, sx=Ld + 2 * nper * a, sy=gy)
    print("saved", f"mode_{tag}.npz")


# FWM coupling of a cavity:  g ~ |int_core (Ep.Ep)(Es.Ei)* dA| / (Np sqrt(Ns Ni)),
# N = int eps |E|^2 dA  (Kerr only in the SiN core).
# Ring reference (travelling mode u(y) e^{i beta x} over L_ring):
#     g_ring = int_core u^4 dy / (int eps u^2 dy)^2 / L_ring
# u(y) is taken from the pump profile at the defect center.   F = g / g_ring.
def check_overlap(args):
    tags = args.tags
    M = {t: np.load(f"mode_{t}.npz") for t in tags}
    out = {}
    for i in range(1, len(tags) - 1):
        ti, tp, ts = tags[i - 1], tags[i], tags[i + 1]
        d = M[tp]; eps = d["eps"]; core = eps > (N_2D**2 + N_OX**2) / 2
        nx, ny = d["Ey"].shape; dA = (float(d["sx"]) / nx) * (float(d["sy"]) / ny); dy = float(d["sy"]) / ny
        Ex = {t: M[t]["Ex"] for t in (ti, tp, ts)}; Ey = {t: M[t]["Ey"] for t in (ti, tp, ts)}
        N = lambda t, vec=True: np.sum(eps * ((abs(Ex[t])**2 if vec else 0) + abs(Ey[t])**2)) * dA
        g_vec = abs(np.sum(core * (Ex[tp]**2 + Ey[tp]**2) * np.conj(Ex[ts]*Ex[ti] + Ey[ts]*Ey[ti])) * dA) \
                / (N(tp) * np.sqrt(N(ts) * N(ti)))
        g_ey = abs(np.sum(core * Ey[tp]**2 * np.conj(Ey[ts] * Ey[ti])) * dA) \
               / (N(tp, False) * np.sqrt(N(ts, False) * N(ti, False)))
        c, w = nx // 2, nx // 20
        col = c - w + np.argmax(np.sum(abs(Ey[tp][c-w:c+w])**2, axis=1))
        u = abs(Ey[tp][col])
        g_ring = np.sum(core[col] * u**4) * dy / (np.sum(eps[col] * u**2) * dy)**2 / L_RING
        out[f"pump {tp} (idler {ti}, signal {ts})"] = dict(F_vector=g_vec / g_ring, F_Ey_only=g_ey / g_ring,
                                                          f=[float(M[t]["fm"]) for t in (ti, tp, ts)])
    print(json.dumps(out, indent=1)); save(out, "overlap.json")


# ----------------------------------------------------------------------------
# PART 4. Summary: tables + net estimate (no MEEP needed)
# ----------------------------------------------------------------------------
def localized_modes(cav, qmin=1e4):
    """Merge the two Harminv points, keep high-Q modes, sort by frequency."""
    seen = {}
    for pt in cav["modes"]:
        for m in pt:
            if m["Q"] > qmin:
                key = round(m["f"], 5)
                seen.setdefault(key, []).append(m["Q"])
    return [(f, min(q), max(q)) for f, q in sorted(seen.items())]


def check_summary(args):
    print("\n=== Check 1: band gap (MPB 3D) ===")
    if os.path.exists("check1.json"):
        d = json.load(open("check1.json"))
        for k, f in d.items():
            delta, a, res = k.split("|")
            print(f"delta {float(delta)*1e3:4.0f} nm  a {a}  {res:>10}:  gap {(f[1]-f[0])*C_THZ:6.3f} THz,"
                  f"  mid-gap {LAM0*F0/(0.5*(f[0]+f[1]))*1e3:7.1f} nm,  next band +{(f[2]-f[1])*C_THZ:5.1f} THz")
    else:
        print("check1.json missing -> run:  gap")

    print("\n=== Check 2: defect modes and F (coupled-mode theory) ===")
    cm = json.load(open("cmt.json")) if os.path.exists("cmt.json") else None
    if cm:
        for k, v in cm.items():
            print(f"{k:9}: L_defect {v['defect_len_um']:.1f} um, decay pump {v['decay_len_pump_um']:.1f} um /"
                  f" signal {v['decay_len_signal_um']:.1f} um,  F = {v['F_vs_ring']:.2f}"
                  f"  (length-only {v['F_length_only']:.1f})")
    else:
        print("cmt.json missing -> run:  cmt")

    F_src = "CMT"
    F = np.mean([v["F_vs_ring"] for v in cm.values()]) if cm else 3.0
    if os.path.exists("overlap.json"):
        ov = json.load(open("overlap.json"))
        print("\n=== F from the simulated fields (MEEP, tapered) ===")
        for k, v in ov.items():
            print(f"{k}:  F = {v['F_Ey_only']:.2f} (Ey only), {v['F_vector']:.2f} (vector)")
        F = np.mean([0.5 * (v["F_Ey_only"] + v["F_vector"]) for v in ov.values()]); F_src = "fields"
    print(f"\n=== Check 3: MEEP defect cavity (2D), Q_int = {Q_INT:.0e}, linewidth ~{LINEWIDTH_GHZ:.2f} GHz,"
          f" F = {F:.1f} ({F_src}) ===")
    print("Required pump vs the plain ring:  P/P_ring = r_p * sqrt(r_s r_i) / F,   r = 1 + Q_int/Q_rad")
    print("(range: worst / best Q_rad of the two probes; valid only if the mismatch is tuned to ~0)")
    for tag in ("abrupt", "taper"):
        fn = f"cav_{tag}.json"
        if not os.path.exists(fn):
            print(f"{fn} missing -> run:  cavity --tag {tag}"); continue
        ms = localized_modes(json.load(open(fn)))
        print(f"\n[{tag}]  modes:")
        for f, q1, q2 in ms:
            print(f"   f {f:.5f}   lambda {1e3/f:7.2f} nm   Q_rad {q1:.1e} - {q2:.1e}")
        for i in range(1, len(ms) - 1):          # every consecutive triplet: idler, pump, signal
            (fi, qi1, qi2), (fp, qp1, qp2), (fs, qs1, qs2) = ms[i-1], ms[i], ms[i+1]
            if fs - fp < 1e-4 or fp - fi < 1e-4:
                continue                         # duplicate Harminv entry
            mism = (fs + fi - 2 * fp) * C_THZ * 1e3
            r = lambda q: 1 + Q_INT / q
            worst = r(qp1) * np.sqrt(r(qs1) * r(qi1)) / F
            best = r(qp2) * np.sqrt(r(qs2) * r(qi2)) / F
            print(f"   pump {1e3/fp:.1f} nm: spacing {(fs-fp)*C_THZ:.2f}/{(fp-fi)*C_THZ:.2f} THz,"
                  f"  mismatch {mism:+7.1f} GHz = {mism/LINEWIDTH_GHZ:+5.0f} linewidths,"
                  f"  P/P_ring = {best:.2f} - {worst:.2f}")
    qneed = Q_INT / (np.sqrt(F / 2) - 1) if F > 2 else float("inf")
    print(f"\nFor P/P_ring = 0.5 with F = {F:.1f}, every mode needs r <= {np.sqrt(F/2):.2f}, i.e. Q_rad >= {qneed:.1e}")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))   # JSON results live next to this script
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gap");  g.add_argument("--deltas", type=float, nargs="+", default=[0.1, 0.2])
    sub.add_parser("cmt")
    c = sub.add_parser("cavity")
    c.add_argument("--delta", type=float, default=0.2, help="corrugation per wall (um)")
    c.add_argument("--Ld", type=float, default=44.4, help="defect length (um)")
    c.add_argument("--Lg", type=float, default=80.0, help="grating length on each side (um)")
    c.add_argument("--ntap", type=int, default=0, help="taper length in periods (0 = abrupt)")
    c.add_argument("--res", type=int, default=25, help="MEEP resolution (px/um)")
    c.add_argument("--T", type=float, default=3000.0, help="run time after the source")
    c.add_argument("--a", type=float, default=0.43987, help="period (um), tuned for 2D mid-gap at 1550")
    c.add_argument("--tag", default="run")
    fd = sub.add_parser("field"); fd.add_argument("fm", type=float); fd.add_argument("tag")
    ov = sub.add_parser("overlap"); ov.add_argument("tags", nargs="+")
    sub.add_parser("summary")
    args = p.parse_args()
    {"gap": check_gap, "cmt": check_cmt, "cavity": check_cavity, "field": check_field,
     "overlap": check_overlap, "summary": check_summary}[args.cmd](args)


if __name__ == "__main__":
    main()
