import numpy as np

def close(x, y, prec=1e-12):
    return np.abs(x - y) < prec

def close_zero(x, prec=1e-12):
    return np.abs(x) < prec

def _round_middle_down(arr):
    arrtmp = np.around(arr)
    mv = (arr - np.floor(arr))==.5
    arrtmp[mv] = np.floor(arr[mv])
    return arrtmp

def periodic_dist(coord1, coord2, torus):
    diff = coord2 - coord1
    if not np.all(close_zero(torus)):
        # remove zero columns
        reduced_torus = torus[~np.all(close_zero(torus), axis=1)].T


        # Compute projection of diff on reduced torus
        # by solving A^TAy=A^Tx ()
        y = np.linalg.solve(np.dot(reduced_torus.T, reduced_torus), 
                            np.dot(reduced_torus.T, diff))

        # Compute ortho
        ortho = diff - np.dot(reduced_torus, y)

        # Pull coordinate back to unit torus
        y -= _round_middle_down(y) 
        proj = np.dot(reduced_torus, y)

        dist = np.linalg.norm(ortho + proj)
    else:
        dist = np.linalg.norm(diff)
    return dist

def get_nb_bonds(nb, coordinates, torus):
    
    # Compute distances
    latticedist = np.array([])
    for c1 in coordinates:
        for c2 in coordinates:
            latticedist = np.append(latticedist,
                                    periodic_dist(c1, c2, torus))
    latticedist = np.unique(np.round(latticedist, decimals=8))
    bonds=[]
    if nb == 0:
        for site, _ in enumerate(coordinates):
            bonds.append(site)
    else:
        for site1, c1 in enumerate(coordinates):
            for site2, c2 in enumerate(coordinates):   
                if np.abs(periodic_dist(c1, c2, torus) - latticedist[nb]) < 1e-8:
                    if site1 < site2:
                        bonds.append([site1, site2])
    return bonds
