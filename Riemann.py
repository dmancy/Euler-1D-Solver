import numpy as np
import math
import matplotlib.pyplot as plt

from Solution import Solution
import Find
from State import State

class Riemann:
    """Definition of a 1D Riemann Problem:
    - W_left       : Left state
    - W_right      : Right state
    - W_left_star  : Left star region
    - W_right_star : Right star region
    - gamma        : Heat capacity ratio
    - P_star       : Pressure in region Left star and Right star
    - U_star       : Velocity in region Left star and Right star"""
    

    def __init__(self, rho_l, u_l, p_l, rho_r, u_r, p_r, gamma): 
            """Constructor of a Riemann Problem"""
            self.W_left = State("Left", gamma, rho_l, u_l, p_l)
            self.W_right = State("Right", gamma, rho_r, u_r, p_r)
            self.gamma = gamma
            self.P_star = Find.Find_P(self.W_left, self.W_right, gamma)
            self.U_star = Find.Find_U(self.W_left, self.W_right, self.P_star, gamma)

            if (self.P_star < self.W_right.pressure):
                #Rarefaction wave
                rho_star_r = self.W_right.rho * (self.P_star/self.W_right.pressure)**(1/gamma)
                self.W_right_star = State("Right star", gamma, rho_star_r, self.U_star, self.P_star)
            else:
                #Shock wave
                rho_star_r = self.W_right.rho * (self.P_star/self.W_right.pressure + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1) * self.P_star/self.W_right.pressure + 1)
                self.W_right_star =  State("Right star", gamma, rho_star_r, self.U_star, self.P_star)


            if (self.P_star < self.W_left.pressure):
                #Rarefaction wave
                rho_star_l = self.W_left.rho * (self.P_star/self.W_left.pressure)**(1/gamma)
                self.W_left_star =  State("Left star", gamma, rho_star_l, self.U_star, self.P_star)
            else:
                #Shock wave
                rho_star_l = self.W_left.rho * (self.P_star/self.W_left.pressure + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1) * self.P_star/self.W_left.pressure + 1)
                self.W_left_star =  State("Left star", gamma, rho_star_l, self.U_star, self.P_star)



    def eval_sampling_point(self, sampling_point):
            """Compute the state at a given sampling point x/t"""
            return Solution(self, sampling_point)

    def plot_time(self, X, x0, t):
            """Plot the repartition of density, velocity and pressure at a given time t, on a mesh X, where the origin of the discontinuity is in x0. Save raw data in a file."""

            Pressure = [Solution(self, (x-x0)/t).pressure for x in X]
            Velocity = [Solution(self, (x-x0)/t).velocity for x in X]
            Density = [Solution(self, (x-x0)/t).rho for x in X]

            fig, axs = plt.subplots(3, sharex=True)
            fig.suptitle("Solution of the Riemann problem\nat t = {}s".format(t))
            axs[0].plot(X, Density)
            axs[1].plot(X, Velocity)
            axs[2].plot(X, Pressure)

            axs[0].grid()
            axs[0].set(ylabel = "Density")
            axs[1].grid()
            axs[1].set(ylabel = "Velocity")
            axs[2].grid()
            axs[2].set(ylabel = "Pressure")

            plt.xlabel("Location x")
