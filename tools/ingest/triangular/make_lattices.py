#!/usr/bin/env python3
"""Generate the 13 triangular.aniso lattices in the library's TOML format.

Geometry comes from the ORIGINAL metts generator (latt.py / create_triangular.py),
which reproduces all 37 surviving triangular lattices in hubbard.finitet and the
metts software tree byte-for-byte, and whose energies match the stored H to
0.07% with no free parameters.

Bond rule (recovered from the data, see CLAUDE.md in the data repo):
    T  on a1 = (1, 0)  and  a2 = (0.5, sqrt3/2)
    Tp on a2 - a1 = (-0.5, sqrt3/2)        -- the diagonal
so tp -> 0 is the square lattice and tp = 1 the isotropic triangular one.

NAMING. The cluster calls the two cylinders `op` and `rect.op`, which does not
say what they are. They are the standard YC and XC cylinders, established by
matching bulk bond offsets against the surviving triangular.{XC4,YC3} lattices
in triangular.heisenberg.dynamics, so the library names them YC and XC. YC/XC
already imply the cylinder (open x, periodic y), so the `op` code is dropped.
`oo` (both open) and `pp` (torus) are not cylinders and keep their codes.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_lat import coordinates_and_torus, coordinates_and_torus_xc

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lattices")
os.makedirs(OUT, exist_ok=True)

TP_DIR = (-0.5, 0.866)          # a2 - a1, rounded as in the class labels

# (library name, cluster name, nx, ny, boundary, cylinder kind)
WANTED = [
    ("triangular.aniso.YC.nx.4.ny.4",  "triangular.aniso.op.nx.4.ny.4",   4, 4, "op", "plain"),
    ("triangular.aniso.YC.nx.8.ny.3",  "triangular.aniso.op.nx.8.ny.3",   8, 3, "op", "plain"),
    ("triangular.aniso.YC.nx.8.ny.4",  "triangular.aniso.op.nx.8.ny.4",   8, 4, "op", "plain"),
    ("triangular.aniso.YC.nx.16.ny.3", "triangular.aniso.op.nx.16.ny.3", 16, 3, "op", "plain"),
    ("triangular.aniso.YC.nx.16.ny.4", "triangular.aniso.op.nx.16.ny.4", 16, 4, "op", "plain"),
    ("triangular.aniso.YC.nx.16.ny.6", "triangular.aniso.op.nx.16.ny.6", 16, 6, "op", "plain"),
    ("triangular.aniso.YC.nx.24.ny.3", "triangular.aniso.op.nx.24.ny.3", 24, 3, "op", "plain"),
    ("triangular.aniso.YC.nx.24.ny.4", "triangular.aniso.op.nx.24.ny.4", 24, 4, "op", "plain"),
    ("triangular.aniso.YC.nx.32.ny.3", "triangular.aniso.op.nx.32.ny.3", 32, 3, "op", "plain"),
    ("triangular.aniso.YC.nx.32.ny.4", "triangular.aniso.op.nx.32.ny.4", 32, 4, "op", "plain"),
    ("triangular.aniso.XC.nx.16.ny.4", "triangular.aniso.rect.op.nx.16.ny.4", 16, 4, "op", "rect"),
    ("triangular.aniso.oo.nx.8.ny.2",  "triangular.aniso.oo.nx.8.ny.2",   8, 2, "oo", "plain"),
    ("triangular.aniso.pp.nx.4.ny.4",  "triangular.aniso.pp.nx.4.ny.4",   4, 4, "pp", "plain"),
]


def displacement(coords, torus, i, j):
    d = coords[j] - coords[i]
    if not np.all(np.abs(torus) < 1e-12):
        red = torus[~np.all(np.abs(torus) < 1e-12, axis=1)].T
        y = np.linalg.solve(red.T @ red, red.T @ d)
        ortho = d - red @ y
        y -= np.round(y)
        d = ortho + red @ y
    return d


def build(nx, ny, boundary, kind):
    coords, torus = (coordinates_and_torus_xc(nx, ny, boundary) if kind == "rect"
                     else coordinates_and_torus(nx, ny, boundary))
    n = nx * ny
    T, Tp = [], []
    for i in range(n):
        for j in range(i + 1, n):
            d = displacement(coords, torus, i, j)
            if abs(np.linalg.norm(d) - 1.0) > 1e-6:
                continue
            v = d if (d[1] > 1e-9 or (abs(d[1]) < 1e-9 and d[0] > 0)) else -d
            key = (round(v[0], 4), round(v[1], 4))
            (Tp if abs(key[0] - TP_DIR[0]) < 1e-3 and abs(key[1] - TP_DIR[1]) < 1e-3
             else T).append((i, j))
    return coords, T, Tp


def to_toml(coords, T, Tp):
    out = ["Coordinates = ["]
    for c in coords:
        out.append(f"  [{c[0]:.12g}, {c[1]:.12g}],")
    out += ["]", "", "Interactions = ["]
    for i, j in T:
        out.append(f"  ['T', 'HUBBARDHOP', {i}, {j}],")
    for i, j in Tp:
        out.append(f"  ['Tp', 'HUBBARDHOP', {i}, {j}],")
    out += ["]", ""]
    return "\n".join(out)


if __name__ == "__main__":
    print(f"{'library name':<34} {'cluster name':<38} {'sites':>6} {'T':>5} {'Tp':>5}  check")
    for name, cluster, nx, ny, bnd, kind in WANTED:
        coords, T, Tp = build(nx, ny, bnd, kind)
        used = {s for b in T + Tp for s in b}
        ok = "ok" if len(used) == nx * ny else f"ONLY {len(used)}/{nx*ny} SITES"
        with open(os.path.join(OUT, name + ".toml"), "w") as f:
            f.write(to_toml(coords, T, Tp))
        print(f"{name:<34} {cluster:<38} {nx*ny:>6} {len(T):>5} {len(Tp):>5}  {ok}")
    print(f"\nwritten to {OUT}")
