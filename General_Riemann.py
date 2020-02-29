import numpy as np
import matplotlib.pyplot as plt

from Solution import Solution
from Riemann import Riemann


def General_Riemann_Problem(W_left, W_right, sampling_point): 
    """return the state (rho, u, p) of the Riemann problem at x/t = sampling_point"""

    Riemann_Problem = Riemann(W_left.density, W_left.velocity, W_left.pressure, W_right.density, W_right.velocity, W_right.pressure, W_right.gamma) 

    W_Solution = Solution(Riemann_Problem, sampling_point)

    return (W_Solution)

