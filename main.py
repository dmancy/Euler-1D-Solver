#Python libraries
import matplotlib.pyplot as plt
import numpy as np

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


gamma = 1.4

U_initial = [State("Left", gamma, 1, 0, 2) if grid.cell_position[i] <= 0 else State("Right", gamma, 1, 0, 1) for i in range(N_cells)]

pressure = [U_initial[i].pressure for i in range(N_cells)]

#plt.figure()
#plt.plot(grid.cell_position, pressure)
#plt.show()

U2 = U_initial.copy()

pressure = [U2[i].pressure for i in range(N_cells)]

#plt.plot(grid.cell_position, pressure)
#plt.show()

t0 = 0
t_final = 25


#Exact solution
Riemann_problem = Riemann(1., 0., 2., 1., 0., 1., 1.4)
#Riemann_problem.plot_time(grid.cell_position, 0,  t_final)

Courant_number = .1

Euler = Euler(U_initial, grid, Courant_number, t0, t_final)
#pressure = [Euler.U_final[i].rho for i in range(N_cells)]
#plt.figure()
#plt.plot(grid.cell_position, pressure, "+")

#U_final = Solver(U_initial, grid, Courant_number, t0, t_final)

plot(Riemann_problem, Euler, grid, t_final)



#pressure = [U_final[i].rho for i in range(N_cells)]
#plt.figure()
#plt.plot(grid.cell_position, pressure, "+")


plt.show()
