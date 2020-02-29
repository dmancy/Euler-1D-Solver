import numpy as np

from General_Riemann import General_Riemann_Problem
from Flux import Flux


def Solver(U_init, grid, Courant_Number, t0 t_final):

    U = np.copy(U_init)

    for i_cell in range(len(grid.cell_position)):
        face1 = grid.faces[i_cell]
        face2 = grid.faces[i_cell+1]

        delta_x = face2 - face1

        delta_t = Courant_Number * delta_x / (abs(U[i].velocity) + U[i].c)

        #Resolution of the Riemann Problems

        W_face_1 = General_Riemann_Problem(U[i-1], U[i])
        W_face_2 = General_Riemann_Problem(U[i], U[i+1])

        Flux1 = Flux(W_face_1)
        Flux2 = Flux(W_face_2)

        #Update values in U
        #Create a method in State to update temperature and c









    
