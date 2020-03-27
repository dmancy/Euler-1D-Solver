#Python libraries
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

from Grid import Grid
from State import State
from Solver import Solver
from Riemann import Riemann
from Euler import Euler
from plot import plot


#Number of cells
N_cells = 100

#Grid step
delta_x = 1

#Number of faces
N_faces = N_cells + 1

#Generation of face positions
faces = np.arange(-N_cells//2, N_cells//2+1, delta_x)


#Grid generation
grid = Grid(faces)


#Specific heat ratio
gamma = 1.4

#Initial condition
U_initial = [State("Left", gamma, 1, 0, 2) if grid.cell_position[i] <= 0 else State("Right", gamma, 1, 0, 1) for i in range(N_cells)]

#Time parameters
t0 = 0
t_final = 25


#Exact solution
Riemann_problem = Riemann(1., 0., 2., 1., 0., 1., 1.4)

Courant_number = 0.7

solution = Euler(U_initial, grid, Courant_number, t0, t_final, "Steger_Warming")

plot(Riemann_problem, solution, grid, t_final, Courant_number)

#Plot Figures
plt.show()
