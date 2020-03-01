import numpy as np

from General_Riemann import General_Riemann_Problem
from Flux import Flux


def Solver(U_init, grid, Courant_Number, t0, t_final):

    U = U_init.copy()
    U_new = U.copy()


    t = t0

    #Uniform Grid
    delta_x = grid.cell_length[0]

    while t <= t_final:
        #Find max(|U| + c)
        max_eigen = 0
        for i_cell in range(len(grid.cell_position)):
            if max_eigen < abs(U[i_cell].velocity) + U[i_cell].c:
                umax = abs(U[i_cell].velocity) + U[i_cell].c

        delta_t = Courant_Number * delta_x / umax

        for i_cell in range(len(grid.cell_position)):
            face1 = grid.faces[i_cell]
            face2 = grid.faces[i_cell+1]

            #Resolution of the Riemann Problems

            W_face_1 = General_Riemann_Problem(U[i_cell-1], U[i_cell], 0)
            W_face_2 = General_Riemann_Problem(U[i_cell], U[(i_cell+1)%len(grid.cell_position)], 0)

            Flux1 = Flux(W_face_1)
            Flux2 = Flux(W_face_2)

            U_new[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)

        U = U_new.copy()

        if (t == t_final):
            break

        t += delta_t

        if ( t > t_final):
            delta_t -= t + delta_t - t_final

        return U
       











    
