import numpy as np

from General_Riemann import General_Riemann_Problem
from Flux import Flux
from Flux_Splitting import Flux_Steger_Warming
from Flux_Van_Leer import Flux_Van_Leer
from Flux_Zha_Bilgen import Flux_Zha_Bilgen
from Flux_Roe import Flux_Roe
import Flux_Richtmeyer

class Euler:

    def __init__(self, U_init, grid, Courant_number, t0, t_final):
        self.U_init = U_init
        self.grid = grid
        self.Courant_Number = Courant_number
        self.t0 = t0
        self.t_final = t_final

        self.U_final = None

        self.Solver_Richtmeyer()


    def Solver_Godunov(self):

        U = self.U_init.copy()

        t = self.t0

        delta_x = min(self.grid.cell_length)

        W_faces = [0] * len(self.grid.faces)

        while t < self.t_final:

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

            t += delta_t
            if ( t > self.t_final):
                delta_t -= t  - self.t_final

            for i_cell in range(len(self.grid.cell_position)):

                W_face_1 = W_faces[i_cell]
                W_face_2 = W_faces[i_cell + 1]
                
                Flux1 = Flux(W_face_1)
                Flux2 = Flux(W_face_2)

                U[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)



            self.U_final = U


    def Solver_Steger_Warming(self):

        U = self.U_init.copy()

        t = self.t0

        delta_x = min(self.grid.cell_length)

        F_faces = [0] * len(self.grid.faces)

        while t < self.t_final:

            #Find max(|U| + c)
            max_eigen = 0
            for i_cell in range(len(self.grid.cell_position)):
                if max_eigen < abs(U[i_cell].velocity) + U[i_cell].c:
                    umax = abs(U[i_cell].velocity) + U[i_cell].c



            delta_t = self.Courant_Number * delta_x / umax

            t += delta_t
            if ( t > self.t_final):
                delta_t -= t  - self.t_final

            #Compute the fluxes at the interfaces
            for i_face in range(len(self.grid.faces)):
                if i_face == 0:
                    Flux = Flux_Steger_Warming(U[i_face],U[i_face])
                elif i_face == len(self.grid.faces)-1:
                    Flux = Flux_Steger_Warming(U[i_face-1],U[i_face-1])

                else:
                    Flux = Flux_Steger_Warming(U[i_face-1], U[i_face])

                F_faces[i_face] = Flux


            #Update the state at each cell
            for i_cell in range(len(self.grid.cell_position)):
                Flux1 = F_faces[i_cell]
                Flux2 = F_faces[i_cell+1]

                U[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)


            self.U_final = U



    def Solver_Van_Leer(self):

        U = self.U_init.copy()

        t = self.t0

        delta_x = min(self.grid.cell_length)

        F_faces = [0] * len(self.grid.faces)


        while t < self.t_final:

            #Find max(|U| + c)
            max_eigen = 0
            for i_cell in range(len(self.grid.cell_position)):
                if max_eigen < abs(U[i_cell].velocity) + U[i_cell].c:
                    umax = abs(U[i_cell].velocity) + U[i_cell].c



            delta_t = self.Courant_Number * delta_x / umax

            t += delta_t
            if ( t > self.t_final):
                delta_t -= t  - self.t_final

            #Compute the fluxes at the interfaces
            for i_face in range(len(self.grid.faces)):
                if i_face == 0:
                    Flux = Flux_Van_Leer(U[i_face],U[i_face])
                elif i_face == len(self.grid.faces)-1:
                    Flux = Flux_Van_Leer(U[i_face-1],U[i_face-1])

                else:
                    Flux = Flux_Van_Leer(U[i_face-1], U[i_face])

                F_faces[i_face] = Flux

            #Update the state at each cell
            for i_cell in range(len(self.grid.cell_position)):
                Flux1 = F_faces[i_cell]
                Flux2 = F_faces[i_cell+1]

                U[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)

            self.U_final = U


    def Solver_Zha_Bilgen(self):

        U = self.U_init.copy()

        t = self.t0

        delta_x = min(self.grid.cell_length)

        F_faces = [0] * len(self.grid.faces)

        while t < self.t_final:

            #Find max(|U| + c)
            max_eigen = 0
            for i_cell in range(len(self.grid.cell_position)):
                if max_eigen < abs(U[i_cell].velocity) + U[i_cell].c:
                    umax = abs(U[i_cell].velocity) + U[i_cell].c



            delta_t = self.Courant_Number * delta_x / umax

            t += delta_t
            if ( t > self.t_final):
                delta_t -= t  - self.t_final

            #Compute the fluxes at the interfaces
            for i_face in range(len(self.grid.faces)):
                if i_face == 0:
                    Flux = Flux_Zha_Bilgen(U[i_face],U[i_face])
                elif i_face == len(self.grid.faces)-1:
                    Flux = Flux_Zha_Bilgen(U[i_face-1],U[i_face-1])

                else:
                    Flux = Flux_Zha_Bilgen(U[i_face-1], U[i_face])

                F_faces[i_face] = Flux

            #Update the state at each cell
            for i_cell in range(len(self.grid.cell_position)):
                Flux1 = F_faces[i_cell]
                Flux2 = F_faces[i_cell+1]

                U[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)

            self.U_final = U


    def Solver_Richtmeyer(self):

        U = self.U_init.copy()

        t = self.t0

        delta_x = min(self.grid.cell_length)

        F_faces = [0] * len(self.grid.faces)

        while t < self.t_final:

            #Find max(|U| + c)
            max_eigen = 0
            for i_cell in range(len(self.grid.cell_position)):
                if max_eigen < abs(U[i_cell].velocity) + U[i_cell].c:
                    umax = abs(U[i_cell].velocity) + U[i_cell].c



            delta_t = self.Courant_Number * delta_x / umax

            t += delta_t
            if ( t > self.t_final):
                delta_t -= t  - self.t_final

            #Compute the Fluxes
            for i_face in range(len(self.grid.faces)):
                if (0 < i_face < len(self.grid.faces)-1): 
                    F_faces[i_face] = Flux(Flux_Richtmeyer.Predictor(U[i_face-1], U[i_face], delta_t, delta_x))
                elif (i_face == 0):
                    F_faces[i_face] = Flux(Flux_Richtmeyer.Predictor(U[i_face], U[i_face], delta_t, delta_x))
                else:
                    F_faces[i_face] = Flux(Flux_Richtmeyer.Predictor(U[i_face-1], U[i_face-1], delta_t, delta_x))



            #Update the state at each cell
            for i_cell in range(len(self.grid.cell_position)):
                Flux1 = F_faces[i_cell]
                Flux2 = F_faces[i_cell+1]

                U[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)


            self.U_final = U

    def Solver_Roe(self):

        U = self.U_init.copy()

        t = self.t0

        delta_x = min(self.grid.cell_length)

        F_faces = [0] * len(self.grid.faces)

        while t < self.t_final:

            #Find max(|U| + c)
            max_eigen = 0
            for i_cell in range(len(self.grid.cell_position)):
                if max_eigen < abs(U[i_cell].velocity) + U[i_cell].c:
                    umax = abs(U[i_cell].velocity) + U[i_cell].c



            delta_t = self.Courant_Number * delta_x / umax

            t += delta_t
            if ( t > self.t_final):
                delta_t -= t  - self.t_final

            #Compute the fluxes at the interfaces
            for i_face in range(len(self.grid.faces)):
                if i_face == 0:
                    Flux = Flux_Roe(U[i_face],U[i_face])
                elif i_face == len(self.grid.faces)-1:
                    Flux = Flux_Roe(U[i_face-1],U[i_face-1])

                else:
                    Flux = Flux_Roe(U[i_face-1], U[i_face])

                F_faces[i_face] = Flux

            #Update the state at each cell
            for i_cell in range(len(self.grid.cell_position)):
                Flux1 = F_faces[i_cell]
                Flux2 = F_faces[i_cell+1]

                U[i_cell].solver_step(delta_t, delta_x, Flux1, Flux2)


            self.U_final = U
