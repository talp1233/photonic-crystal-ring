# Vector perturbation theory for the Bragg splitting of a sidewall corrugation.
# Standing waves at the band edge: transverse E ~ cos(bz), longitudinal Ez ~ sin(bz),
# so Ez enters the splitting with the OPPOSITE sign to the transverse parts.
import meep as mp
from meep import mpb
import numpy as np, json
f0 = 1 / 1.55; n_sin, n_ox = 1.996, 1.444; W0, T0 = 1.0, 0.8
sx = sy = 4.0; res = 64
geo = [mp.Block(size=mp.Vector3(W0, T0, mp.inf), material=mp.Medium(index=n_sin))]
ms = mpb.ModeSolver(geometry_lattice=mp.Lattice(size=mp.Vector3(sx, sy, 0)), geometry=geo,
                    default_material=mp.Medium(index=n_ox), resolution=res, num_bands=1,
                    k_points=[mp.Vector3(0, 0, 1.14427)])
ms.run()
f = ms.get_freqs()[0]
E = np.array(ms.get_efield(1, bloch_phase=False))[:, :, 0, :]
eps = np.array(ms.get_epsilon())
nx, ny = eps.shape; dx = sx / nx; dy = sy / ny
xs = (np.arange(nx) + 0.5) * dx - sx / 2; ys = (np.arange(ny) + 0.5) * dy - sy / 2
D = eps[..., None] * E
den = np.sum(eps * np.sum(np.abs(E) ** 2, -1)) * dx * dy
de = n_sin ** 2 - n_ox ** 2; dinv = 1 / n_sin ** 2 - 1 / n_ox ** 2
rows = np.where(np.abs(ys) < T0 / 2)[0]
def wall_integral(sign_z):
    tot = 0
    for xw in (-W0 / 2, W0 / 2):
        i_in = np.argmin(np.abs(xs - (xw - np.sign(xw) * dx)))   # just inside the core
        i_out = np.argmin(np.abs(xs - (xw + np.sign(xw) * dx)))  # just outside
        for j in rows:
            Ey = 0.5 * (abs(E[i_in, j, 1]) ** 2 + abs(E[i_out, j, 1]) ** 2)  # continuous
            Ez = 0.5 * (abs(E[i_in, j, 2]) ** 2 + abs(E[i_out, j, 2]) ** 2)
            Dx = 0.5 * (abs(D[i_in, j, 0]) ** 2 + abs(D[i_out, j, 0]) ** 2)
            tot += (de * (Ey + sign_z * Ez) - dinv * Dx) * dy
    return tot
uniform = 0.5 * wall_integral(+1) / den        # |d(omega)/omega| per um of each-wall shift
split = 0.5 * abs(wall_integral(-1)) / den     # splitting per um of sinusoidal amplitude
nu0 = 3e8 / 1.55e-6
out = dict(freq=float(f), uniform_rel_per_um=float(uniform),
           uniform_check_from_dfdW=float(2 * 0.05375979991349422 / f0),
           split_GHz_per_nm=float(split * 1e-3 * nu0 / 1e9),
           Ez_share_of_wall_integrand=float(abs(wall_integral(-1)) / wall_integral(+1)))
print(json.dumps(out, indent=1))
