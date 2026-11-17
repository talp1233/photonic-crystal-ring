# Direct 3D check: corrugated waveguide unit cell, gap at the Bragg point
import meep as mp
from meep import mpb
import numpy as np, sys, json
n_sin, n_ox = 1.996, 1.444
W0, T0 = 1.0, 0.8
a = 0.43696                      # Bragg period from the cross-section run (um)
res = int(sys.argv[2]) if len(sys.argv) > 2 else 40
delta = float(sys.argv[1])       # sinusoidal amplitude per wall (um), both walls
ns = int(sys.argv[4]) if len(sys.argv)>4 else 32
geo = []
for i in range(ns):
    x = (i + 0.5) / ns * a - a / 2
    w = W0 + 2 * delta * np.cos(2 * np.pi * x / a)
    geo.append(mp.Block(center=mp.Vector3(x, 0, 0), size=mp.Vector3(a / ns, w, T0),
                        material=mp.Medium(index=n_sin)))
ms = mpb.ModeSolver(geometry_lattice=mp.Lattice(size=mp.Vector3(a, 3.0, 3.0)),
                    geometry=geo, default_material=mp.Medium(index=n_ox),
                    resolution=mp.Vector3(int(sys.argv[3]) if len(sys.argv)>3 else res, res, res), num_bands=3, k_points=[mp.Vector3(0.5, 0, 0)])
ms.run_zeven()
f = list(ms.get_freqs())
print("RESULT", delta, res, json.dumps([float(v) for v in f]))
