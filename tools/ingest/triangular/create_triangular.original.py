#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from latt import *

ny=4
nxs=[4, 8, 16]
xperiodic = False
yperiodic = True
plot=True
numbering = "yfirst" #either "snake" or "yfirst"

for nx in nxs:

    coordinates = np.zeros((nx*ny, 2))
    cell = np.array([[nx, 0], [0.5*ny, np.sqrt(3)/2*ny]])
    unit_cell = np.array([[1, 0],[0.5, np.sqrt(3)/2]])
    basis = np.array([[0, 0]])

    if xperiodic and yperiodic:
        torus = cell
        periodic_string = "pp"
    elif xperiodic:
        torus = np.array([[nx, 0],[0, 0]])
        periodic_string = "po"
    elif yperiodic:
        torus = np.array([[0, 0], [0.5*ny, np.sqrt(3)/2*ny]])
        periodic_string = "op"
    else:
        torus = np.array([[0, 0],[0, 0]])
        periodic_string = "oo"

    filename = "triangular.nx.{}.ny.{}.{}.hubbard.hop.{}.lat".format(nx, ny, numbering, periodic_string)

    # Create snake geometry
    idx = 0
    for a in range(nx+ny):
        for b in range(ny):
            if numbering == "snake":
                x = a - b
                y = b
            elif numbering == "yfirst":
                x = a
                y = b

            if (0 <= x) and (x < nx):
                for bc in basis:
                    coordinates[idx, :] = x*unit_cell[0, :] + y*unit_cell[1, :] + bc
                    idx += 1


    nn_hoppings = get_nb_bonds(1, coordinates, torus)
    nnn_hoppings = get_nb_bonds(2, coordinates, torus)

    for hop in nn_hoppings:
        print("nn", hop)
    for hop in nnn_hoppings:
        print("nnn", hop)


    if plot:
        fig, ax = plt.subplots(1)
        ax.plot(coordinates[:,0], coordinates[:,1], "o")
        for i,c in enumerate(coordinates):
            ax.text(c[0], c[1], i)

        for hop in nn_hoppings:
            ax.plot([coordinates[hop[0], 0], coordinates[hop[1], 0]],
                    [coordinates[hop[0], 1], coordinates[hop[1], 1]], 
                    "-", color="purple")
        for hop in nnn_hoppings:
            ax.plot([coordinates[hop[0], 0], coordinates[hop[1], 0]],
                    [coordinates[hop[0], 1], coordinates[hop[1], 1]], 
                    "-", color="firebrick")
        plt.show()

    with open(filename, 'w') as fl:
        fl.write("[Interactions]\n")
        for hop in nn_hoppings:
            fl.write("HUBBARDHOP T {} {}\n".format(hop[0], hop[1]))
        for hop in nnn_hoppings:
            fl.write("HUBBARDHOP Tp {} {}\n".format(hop[0], hop[1]))
    print("Written lattice to file ", filename)

