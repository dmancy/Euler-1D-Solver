import numpy as np
import matplotlib.pyplot as plt

def function_F(W, gamma, P_eval):
    #Evaluation of the function f inplied in the search of P_star and U_star
    if (P_eval > W.pressure):
        return (P_eval - W.pressure) * (2/((gamma+1) * W.rho) / (P_eval + (gamma-1)/(gamma+1) * W.pressure))**.5
    else:
        return (2*W.c/(gamma-1)) * ((P_eval/W.pressure)**((gamma-1)/(2*gamma)) - 1)


def derivative_F(W, gamma, P_eval):
    #Evaluation of the derivative of f used in the Newton method
    if (P_eval > W.pressure):
        return (2/((gamma+1) * W.rho) / (P_eval + (gamma-1)/(gamma+1) * W.pressure))**.5 * (1 - (P_eval - W.pressure)/(2*(P_eval + (gamma-1)/(gamma+1) * W.pressure)))
    else:
        return  1/(W.rho * W.c) * (P_eval/W.pressure)**(-(gamma+1)/(2*gamma))


def Find_P(W_left, W_right, gamma):
    #Find P_star around the contact surface

    #Pressure from the acoustic approximation (Used as the first guess for the Newton method)
    P0 = max((W_left.pressure * W_right.rho * W_right.c + W_right.pressure * W_left.rho * W_left.c + (W_left.velocity - W_right.velocity) * W_right.rho * W_right.c * W_left.rho * W_left.c)/(W_left.rho * W_left.c + W_right.rho * W_right.c),0.000001)

    epsilon = 1
    count = 0

    #Newton method
    while (epsilon > 10**-6 and count < 1000):
        P1 = P0 - (W_right.velocity - W_left.velocity + function_F(W_left, gamma, P0) + function_F(W_right, gamma, P0))/(derivative_F(W_left, gamma, P0) + derivative_F(W_right, gamma, P0))
        epsilon = abs(P1 - P0)/P0
        P0 = P1
        count+= 1
    return P1


def Find_U(W_left, W_right, P_star, gamma):
    #Compute U around the contact surface
    return .5*(W_left.velocity + W_right.velocity) + .5*(function_F(W_right, gamma, P_star) - function_F(W_left, gamma, P_star))

