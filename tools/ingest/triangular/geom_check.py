#!/usr/bin/env python3
"""At tp=1, is `aniso` the nearest-neighbour triangular lattice or the nn+nnn one?

Both aniso readings degenerate to all-three-nn at tp=1, so this does not test
the T/Tp split -- it tests the GEOMETRY, i.e. whether the aniso lattices are
the same bond set as t1t2 (which also puts Tp on the second shell). With U=0
the prediction is parameter-free:

    nn only  :  H = -2 * sum_{shell 1} 2 Re<c+c>
    nn + nnn :  H = -2 * (sum_{shell 1} + sum_{shell 2}) 2 Re<c+c>
"""
import os, re, glob, sys, h5py, numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_lat import coordinates_and_torus, bonds_t1t2

V = ("/data/condmat/awietek/flatiron/ceph/Research/Projects/hubbard.triangular.metts.v2/"
     "metts")

coords, torus = coordinates_and_torus(16, 4, "op")
nn, nnn = bonds_t1t2(coords, torus)
nn = np.array(nn)
nnn = np.array(nnn)
print(f"shell 1: {len(nn)} bonds    shell 2: {len(nnn)} bonds")

files = sorted(glob.glob(
    V + "/**/triangular.aniso.op.nx.16.ny.4*/t.1.tp.1.U.0.holes.0/*/*.h5",
    recursive=True))[:40]
print(f"tp=1, U=0 dumps: {len(files)}")

S1, S2, Hs = [], [], []
for p in files:
    try:
        with h5py.File(p, "r") as h:
            M = np.array(h["CdagupCup"])
            H = np.array(h["H"]).ravel()
    except Exception:
        continue
    ns = min(len(H), M.shape[0])
    if ns == 0:
        continue
    S1.append(2.0 * M[:ns, nn[:, 0], nn[:, 1]].real.sum(axis=1))
    S2.append(2.0 * M[:ns, nnn[:, 0], nnn[:, 1]].real.sum(axis=1))
    Hs.append(H[:ns])

S1 = np.concatenate(S1); S2 = np.concatenate(S2); Hs = np.concatenate(Hs)
print(f"pooled samples: {len(Hs)}\n")

for name, pred in (("nn only   ", -2.0 * S1),
                   ("nn + nnn  ", -2.0 * (S1 + S2))):
    rms = float(np.sqrt(np.mean((Hs - pred) ** 2)))
    # also the best-fit scale, which should be 1.0 for the right geometry
    lam = float(np.dot(pred, Hs) / np.dot(pred, pred))
    print(f"   {name}  rms={rms:9.4f}   best-fit scale={lam:.4f}  (should be 1.0000)")

print(f"\n   mean H = {Hs.mean():.4f},  mean -2*S1 = {(-2*S1).mean():.4f},"
      f"  mean -2*(S1+S2) = {(-2*(S1+S2)).mean():.4f}")
