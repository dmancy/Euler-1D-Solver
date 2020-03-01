import numpy as np

from General_Riemann import General_Riemann_Problem
from Flux import Flux

class Euler:

    def __init__(self, U_init, grid, Courant_number, t0, t_final):
        self.U_init = U_init
        self.grid = grid
        self.Courant_Number = Courant_number
        self.t0 = t0
        self.t_final = t_final

        self.U_final = None

        self.Solver_Godunov()


    def Solver_Godunov(self):

        U = self.U_init.copy()
        U_new = U.copy()

        t = self.t0

        delta_x = min(self.grid.cell_length)

        W_faces = [0] * len(self.grid.faces)

        while t <= self.t_final:
            #Find max(|U| + c)
            max_eigen = 0
            for i_cell in range(len(self.grid.cell_position)):
                if max_eigen < abs(U[i_cell].velocity) + U[i_cell].c:
                    umax = abs(U[i_cell].velocity) + U[i_cell].c

                #Resolution of the Riemann Problems
                if (i_cell != 0):
                    W_faces[i_cell] = General_Riemann_Problem(U[i_cell-1], U[i_cell], 0)
                else:
                    W_faces[i_cell] = U[i_cell]

                W_faces[-1] = U[-1]

            delta_t = self.Courant_Number * delta_x / umax

            for i_cell in range(len(self.grid.cell_position)):

                W_face_1 = W_faces[i_cell]
                W_face_2 = W_faces[i_cell + 1]
                
                Flux1 = Flux(W_face_1)
                Flux2 = Flux(W_face_2)

                U_new[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)

            U = U_new.copy()

            if (t == self.t_final):
                break

            t += delta_t

            if ( t > self.t_final):
                delta_t -= t + delta_t - self.t_final

            self.U_final = U


