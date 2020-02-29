#Python libraries
import matplotlib.pyplot as plt
import numpy as np

from Grid import Grid
from State import State


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

#plt.plot(grid.cell_position, pressure)
#plt.show()




Courant_number = .6
