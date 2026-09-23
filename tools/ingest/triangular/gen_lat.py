#!/usr/bin/env python3
"""Regenerate triangular lattices using the ORIGINAL metts generator machinery.

`latt.py` (periodic_dist / get_nb_bonds) is imported untouched from the metts
software tree, and the coordinate construction is copied verbatim from
create_triangular.py. Nothing here reimplements the geometry -- that is the
point: the t1t2 lattices it produces are diffed against the files that
survived, so the coordinate, torus and distance-shell logic is proven before
the same machinery is used for the aniso bond rule, which has no surviving
example.

Two bond conventions:
  t1t2  -- T on distance shell 1, Tp on shell 2   (what create_triangular.py does)
  aniso -- T on the square bonds (+x, +y), Tp on ONE diagonal; interpolates
           square (tp=0) to isotropic triangular (tp=1)
"""
import sys, os
import numpy as np

sys.path.insert(0, "/home/awietek/Research/Software/metts/old/runtimes/lattice-files")
from latt import get_nb_bonds, periodic_dist          # noqa: E402  (original code)


def coordinates_and_torus(nx, ny, boundary, numbering="yfirst"):
    """Verbatim from create_triangular.py, with the boundary made a parameter."""
    coordinates = np.zeros((nx * ny, 2))
    cell = np.array([[nx, 0], [0.5 * ny, np.sqrt(3) / 2 * ny]])
    unit_cell = np.array([[1, 0], [0.5, np.sqrt(3) / 2]])
    basis = np.array([[0, 0]])

    xperiodic = boundary in ("pp", "po")
    yperiodic = boundary in ("pp", "op")

    if xperiodic and yperiodic:
        torus = cell
    elif xperiodic:
        torus = np.array([[nx, 0], [0, 0]])
    elif yperiodic:
        torus = np.array([[0, 0], [0.5 * ny, np.sqrt(3) / 2 * ny]])
    else:
        torus = np.array([[0, 0], [0, 0]])

    idx = 0
    for a in range(nx + ny):
        for b in range(ny):
            if numbering == "snake":
                x, y = a - b, b
            else:
                x, y = a, b
            if 0 <= x < nx:
                for bc in basis:
                    coordinates[idx, :] = x * unit_cell[0, :] + y * unit_cell[1, :] + bc
                    idx += 1
    assert idx == nx * ny, f"built {idx} sites, expected {nx*ny}"
    return coordinates, torus


def coordinates_and_torus_xc(nx, ny, boundary):
    """XC cylinder (the `.rect` files), read off the surviving XC4 coordinates.

    A column runs up the y axis zig-zagging between a2 = (1/2, sqrt3/2) and
    a2 - a1 = (-1/2, sqrt3/2), so two sites advance by (0, sqrt3) = 2*a2 - a1;
    columns are spaced by a1 = (1, 0). The circumference is therefore vertical,
    (0, ny*sqrt3/2), which is a lattice vector only for even ny -- which is why
    `.rect` exists at ny = 4 and 6 but never 3.
    """
    if ny % 2 != 0:
        raise ValueError(f"XC needs even ny, got {ny}")
    a1 = np.array([1.0, 0.0])
    a2 = np.array([0.5, np.sqrt(3) / 2])

    coordinates = np.zeros((nx * ny, 2))
    idx = 0
    for x in range(nx):
        p = x * a1
        for y in range(ny):
            coordinates[idx, :] = p
            idx += 1
            p = p + (a2 if y % 2 == 0 else a2 - a1)
    assert idx == nx * ny

    xperiodic = boundary in ("pp", "po")
    yperiodic = boundary in ("pp", "op")
    circ = np.array([0.0, ny * np.sqrt(3) / 2])
    if xperiodic and yperiodic:
        torus = np.array([[nx, 0.0], circ])
    elif xperiodic:
        torus = np.array([[nx, 0.0], [0.0, 0.0]])
    elif yperiodic:
        torus = np.array([[0.0, 0.0], circ])
    else:
        torus = np.zeros((2, 2))
    return coordinates, torus


def bonds_t1t2(coords, torus):
    return get_nb_bonds(1, coords, torus), get_nb_bonds(2, coords, torus)


def bonds_aniso(coords, torus, nx, ny, diagonal):
    """T on the square bonds, Tp on one diagonal.

    Selection is by lattice vector, not by distance: in the isotropic
    embedding all three nearest-neighbour directions are the same length, so a
    distance shell cannot separate them. `diagonal` is "+x-y" or "+x+y".
    """
    def site(x, y):
        return x * ny + y

    xper = torus[0, 0] != 0
    yper = not np.all(np.abs(torus[1]) < 1e-12)

    T, Tp = [], []
    for x in range(nx):
        for y in range(ny):
            s = site(x, y)
            # +y
            if y + 1 < ny:
                T.append([s, site(x, y + 1)])
            elif yper and ny > 2:
                T.append(sorted([s, site(x, 0)]))
            # +x
            if x + 1 < nx:
                T.append([s, site(x + 1, y)])
            elif xper and nx > 2:
                Tp_or_T = sorted([s, site(0, y)])
                T.append(Tp_or_T)
            # diagonal
            dy = -1 if diagonal == "+x-y" else 1
            nxp, nyp = x + 1, y + dy
            if nxp < nx and 0 <= nyp < ny:
                Tp.append(sorted([s, site(nxp, nyp)]))
            elif nxp < nx and yper and ny > 2:
                Tp.append(sorted([s, site(nxp, nyp % ny)]))
    dedup = lambda bs: sorted({tuple(b) for b in bs})
    return [list(b) for b in dedup(T)], [list(b) for b in dedup(Tp)]


def write_lat(path, T, Tp):
    with open(path, "w") as f:
        f.write("[Interactions]\n")
        for b in T:
            f.write("HUBBARDHOP T {} {}\n".format(b[0], b[1]))
        for b in Tp:
            f.write("HUBBARDHOP Tp {} {}\n".format(b[0], b[1]))


FT_DIR = ("/data/condmat/awietek/flatiron/ceph/Research/Projects/hubbard.finitet/"
          "lattice-files/")
METTS_DIR = "/home/awietek/Research/Software/metts/old/runtimes/lattice-files/"
TEST_DIR = "/home/awietek/Research/Software/metts/test/old/hubbard/lattice-files/"

if __name__ == "__main__":
    # --- verification: reproduce the surviving t1t2-style files exactly
    SURV = [
        ("/data/condmat/awietek/flatiron/ceph/Research/Projects/hubbard.finitet/lattice-files/"
         "triangular.nx.{}.ny.{}.yfirst.hubbard.hop.{}.lat"),
        ("/home/awietek/Research/Software/metts/old/runtimes/lattice-files/"
         "triangular.nx.{}.ny.{}.yfirst.hubbard.hop.{}.lat"),
    ]
    print("=== regenerating t1t2 lattices and diffing against surviving files")
    ok = bad = miss = 0
    for nx in (4, 8, 16, 32, 64, 128):
        for ny in (2, 3, 4, 6):
            for bnd in ("op", "oo"):
                src = None
                for pat in SURV:
                    p = pat.format(nx, ny, bnd)
                    if os.path.exists(p):
                        src = p
                        break
                if src is None:
                    miss += 1
                    continue
                coords, torus = coordinates_and_torus(nx, ny, bnd)
                T, Tp = bonds_t1t2(coords, torus)
                got = "[Interactions]\n" + "".join(
                    "HUBBARDHOP T {} {}\n".format(*b) for b in T) + "".join(
                    "HUBBARDHOP Tp {} {}\n".format(*b) for b in Tp)
                want = open(src).read()
                if got == want:
                    ok += 1
                    print(f"   MATCH  nx={nx:<4} ny={ny} {bnd}  (T={len(T)}, Tp={len(Tp)})")
                else:
                    bad += 1
                    gl, wl = got.splitlines(), want.splitlines()
                    print(f"   DIFFER nx={nx:<4} ny={ny} {bnd}  got {len(gl)} lines, file {len(wl)}")
    print(f"\n   matched {ok}, differed {bad}, no surviving file {miss}")

    # --- verification: the XC construction against the surviving .rect files
    print("\n=== regenerating XC (.rect) lattices and diffing against surviving files")
    ok = bad = miss = 0
    for nx in (4, 8, 16):
        for ny in (4, 6):
            src = None
            for d in (FT_DIR, METTS_DIR, TEST_DIR):
                p = d + f"triangular.rect.nx.{nx}.ny.{ny}.yfirst.hubbard.hop.op.lat"
                if os.path.exists(p):
                    src = p
                    break
            if src is None:
                miss += 1
                continue
            coords, torus = coordinates_and_torus_xc(nx, ny, "op")
            T, Tp = bonds_t1t2(coords, torus)
            got = "[Interactions]\n" + "".join(
                "HUBBARDHOP T {} {}\n".format(*b) for b in T) + "".join(
                "HUBBARDHOP Tp {} {}\n".format(*b) for b in Tp)
            want = open(src).read()
            if got == want:
                ok += 1
                print(f"   MATCH  nx={nx:<4} ny={ny} op  (T={len(T)}, Tp={len(Tp)})")
            else:
                bad += 1
                wT = sum(1 for l in want.splitlines() if l.split()[1:2] == ["T"])
                wTp = sum(1 for l in want.splitlines() if l.split()[1:2] == ["Tp"])
                print(f"   DIFFER nx={nx:<4} ny={ny} op  got T={len(T)},Tp={len(Tp)}  "
                      f"file T={wT},Tp={wTp}")
    print(f"\n   matched {ok}, differed {bad}, no surviving file {miss}")
