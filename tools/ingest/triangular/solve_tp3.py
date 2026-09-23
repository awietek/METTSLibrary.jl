#!/usr/bin/env python3
"""Same test with the spin factor FIXED at its physical value of 2.

H = -2 * sum_d t_d * sum_{(i,j) in d} 2 Re<c+_i c_j>     (U = 0)

is a zero-parameter prediction once the assignment is chosen, so the residual
is a clean model-selection statistic. The free-scale version could hide a wrong
assignment by inflating the scale; this cannot.
"""
import os, re, glob, sys, h5py, numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solve_tp import nn_classes

V = ("/data/condmat/awietek/flatiron/ceph/Research/Projects/hubbard.triangular.metts.v2/"
     "metts")

files = sorted(glob.glob(
    V + "/**/triangular.aniso.op.nx.16.ny.4*/t.1.tp.*.U.0.holes.0/*/*.h5",
    recursive=True))
bytp = defaultdict(list)
for p in files:
    m = re.search(r"\.tp\.([0-9.]+)\.U\.0\.", p)
    if m:
        bytp[float(m.group(1))].append(p)

cls = nn_classes(16, 4, "plain", "op")
keys = sorted(cls, key=lambda k: (-len(cls[k]), k))
LABEL = {(0.5, 0.866): "a2  (+y)", (1.0, 0.0): "a1  (+x)", (-0.5, 0.866): "a2-a1 (diag)"}

print("class sizes:", {LABEL.get(k, str(k)): len(cls[k]) for k in keys})
print("\nfixed spin factor 2, no free parameters\n")
print(f"{'tp':>5} {'samples':>8}   " + "".join(f"{LABEL.get(k,str(k)):>16}" for k in keys))

totals = defaultdict(float)
for tp in sorted(bytp):
    if abs(tp - 1.0) < 1e-9:
        continue
    S, Hs = [], []
    for p in bytp[tp][:40]:
        try:
            with h5py.File(p, "r") as h:
                M = np.array(h["CdagupCup"])
                H = np.array(h["H"]).ravel()
        except Exception:
            continue
        ns = min(len(H), M.shape[0])
        if ns == 0:
            continue
        row = np.zeros((ns, len(keys)))
        for d, k in enumerate(keys):
            idx = np.array(cls[k])
            row[:, d] = 2.0 * M[:ns, idx[:, 0], idx[:, 1]].real.sum(axis=1)
        S.append(row)
        Hs.append(H[:ns])
    if not S:
        continue
    S, Hs = np.vstack(S), np.concatenate(Hs)

    cells = []
    for c in range(len(keys)):
        t = np.ones(len(keys))
        t[c] = tp
        pred = -2.0 * (S * t).sum(axis=1)
        rms = float(np.sqrt(np.mean((Hs - pred) ** 2)))
        cells.append(rms)
        totals[keys[c]] += rms
    best = int(np.argmin(cells))
    row = "".join(f"{v:>16.4f}" + ("*" if i == best else " ")[:0] for i, v in enumerate(cells))
    mark = "   <- best: " + LABEL.get(keys[best], str(keys[best]))
    print(f"{tp:>5} {len(Hs):>8}   " + "".join(f"{v:>16.4f}" for v in cells) + mark)

print("\nsummed rms over all tp:")
for k, v in sorted(totals.items(), key=lambda kv: kv[1]):
    print(f"   {LABEL.get(k, str(k)):<16} {v:10.4f}")
