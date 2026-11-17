# Cross-section of the Vaidya SiN waveguide (1.0 um wide x 0.8 um thick, SiO2 clad)
# MPB 2D cell (x = width, y = thickness), propagation along z (out of plane k).
import meep as mp
from meep import mpb
import numpy as np, json

lam0 = 1.55
f0 = 1 / lam0                      # frequency in c/um
n_sin, n_ox = 1.996, 1.444
W0, T0 = 1.0, 0.8
sx, sy, res = 4.0, 4.0, 32

def solver(W, T=T0, kz=None, nb=2):
    geo = [mp.Block(size=mp.Vector3(W, T, mp.inf), material=mp.Medium(index=n_sin))]
    ms = mpb.ModeSolver(geometry_lattice=mp.Lattice(size=mp.Vector3(sx, sy, 0)),
                        geometry=geo, default_material=mp.Medium(index=n_ox),
                        resolution=res, num_bands=nb,
                        k_points=[mp.Vector3(0, 0, kz if kz else 1.5)])
    return ms

# 1) propagation constant at 1550 nm, band 1 (fundamental, highest n_eff)
ms = solver(W0)
kz = ms.find_k(mp.NO_PARITY, f0, 1, 1, mp.Vector3(0, 0, 1), 1e-6, f0 * 1.8, f0 * 1.45, f0 * 2.0)[0]
neff = kz / f0
# group index from find_k at f0*(1 +- 0.002)
k_hi = solver(W0).find_k(mp.NO_PARITY, f0 * 1.002, 1, 1, mp.Vector3(0, 0, 1), 1e-7, kz, f0, f0 * 2.1)[0]
k_lo = solver(W0).find_k(mp.NO_PARITY, f0 * 0.998, 1, 1, mp.Vector3(0, 0, 1), 1e-7, kz, f0 * 0.9, f0 * 2.1)[0]
ng = (k_hi - k_lo) / (f0 * 0.004)

# 2) field of band 1 at kz -> polarisation, A_eff
ms = solver(W0, kz=kz); ms.run()
E = ms.get_efield(1, bloch_phase=False)[:, :, 0, :]       # (nx, ny, 3)
eps = ms.get_epsilon()                                     # (nx, ny)
I = np.sum(np.abs(E) ** 2, axis=-1)
pol = [float(np.sum(np.abs(E[..., i]) ** 2) / np.sum(I)) for i in range(3)]
core = eps > (n_ox ** 2 + n_sin ** 2) / 2
dA = (sx / E.shape[0]) * (sy / E.shape[1])                 # um^2 per pixel
Aeff = (np.sum(I) * dA) ** 2 / (np.sum(I[core] ** 2) * dA)  # Kerr only in the core
n2 = 2.4e-19                                               # m^2/W, assumed SiN value
omega = 2 * np.pi * 3e8 / 1.55e-6
gamma = n2 * omega / (3e8 * Aeff * 1e-12)

# 3) d(omega)/dW at fixed kz (both walls move) -> Bragg splitting per corrugation depth
dW = 0.02
ws = []
for W in (W0 - dW, W0 + dW):
    m2 = solver(W, kz=kz); m2.run(); ws.append(m2.get_freqs()[0])
dfdW = (ws[1] - ws[0]) / (2 * dW)             # (c/um) per um
S = 2 * abs(dfdW) / f0                        # relative shift per um of ONE wall's displacement, both walls
nu0 = 3e8 / 1.55e-6
lw_vaidya = nu0 / 2e5                         # loaded linewidth of the Vaidya ring (Hz)
a = lam0 / (2 * neff)                         # Bragg period (um)
L = 2 * np.pi * 120.0
out = dict(neff=neff, ng=ng, polarisation_xyz=pol, Aeff_um2=Aeff, gamma_perWm=gamma,
           dfdW=dfdW, bragg_period_um=a, periods_in_ring=L / a,
           split_GHz_per_nm_bothwalls=S * 1e-3 * nu0 / 1e9,
           split_GHz_per_nm_one_wall=S / 2 * 1e-3 * nu0 / 1e9,
           vaidya_linewidth_GHz=lw_vaidya / 1e9)
out = {k: (float(v) if not isinstance(v, list) else [float(x) for x in v]) for k, v in out.items()}
out["depth_nm_for_10GHz_bothwalls"] = 10 / out["split_GHz_per_nm_bothwalls"]
out["depth_nm_for_10GHz_one_wall"] = 10 / out["split_GHz_per_nm_one_wall"]
out = {k: (float(v) if not isinstance(v, list) else [float(x) for x in v]) for k, v in out.items()}
print(json.dumps(out, indent=1))
json.dump(out, open("xsec.json", "w"), indent=1)
