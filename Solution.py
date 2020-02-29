import math
import numpy as np
from State import State


def Solution(R, sampling_point):
        #return the state of a given Riemann problem, at a given sampling point x/t

        gamma = R.W_right.gamma

        if (sampling_point < R.U_star):
            #We are on the left from CS

            if (R.P_star < R.W_left.pressure):
                #Left expansion fan
                c_star_l  = R.W_left.c * (R.P_star/R.W_left.pressure)**((gamma-1)/(2*gamma))
                Speed_HL = R.W_left.velocity - R.W_left.c
                Speed_TL = R.U_star - c_star_l

                if (sampling_point < Speed_HL):
                    #We are before the expansion fan
                    return R.W_left

                elif (sampling_point > Speed_TL):
                    #We are after the expansion fan
                    return R.W_left_star

                else:
                    #We are inside the expansion fan
                    rho = R.W_left.rho * (2/(gamma+1) + (gamma-1)/(gamma+1) * 1/R.W_left.c * (R.W_left.velocity - sampling_point))**(2/(gamma-1))
                    u =  2/(gamma+1) * (  R.W_left.c + (gamma-1)/2*R.W_left.velocity + sampling_point)
                    p = R.W_left.pressure * (2/(gamma+1) + (gamma-1)/(gamma+1) * 1/R.W_left.c * (R.W_left.velocity - sampling_point))**(2*gamma/(gamma-1))

                    return State("Expansion Fan", gamma, rho, u, p)

            else:
                #Left Shock wave
                Speed_Shock = R.W_left.velocity - R.W_left.c * ((gamma+1)/(2*gamma) * (R.P_star/R.W_left.pressure) + (gamma-1)/(2*gamma))**.5

                if(sampling_point < Speed_Shock):
                    #W_left
                    return R.W_left

                else:
                    #W_left_star
                    return R.W_left_star

        else:
            #We are on the right from CS

            if (R.P_star < R.W_right.pressure):
                #Right fan
                c_star_r  = R.W_right.c * (R.P_star/R.W_right.pressure)**((gamma-1)/(2*gamma))
                Speed_HR = R.W_right.velocity + R.W_right.c
                Speed_TR = R.U_star + c_star_r

                if (sampling_point > Speed_HR):
                    return R.W_right

                elif (sampling_point < Speed_TR):
                    return R.W_right_star

                else:
                    rho = R.W_right.rho * (2/(gamma+1) - (gamma-1)/(gamma+1) * 1/R.W_right.c * (R.W_right.velocity - sampling_point))**(2/(gamma-1))
                    u =  2/(gamma+1) * (  -R.W_right.c + (gamma-1)/2*R.W_right.velocity + sampling_point)
                    p = R.W_right.pressure * (2/(gamma+1) - (gamma-1)/(gamma+1) * 1/R.W_right.c * (R.W_right.velocity - sampling_point))**(2*gamma/(gamma-1))

                    return State("Expansion Fan", gamma, rho, u, p)

            else:
                #Right Shock Wave
                Speed_Shock = R.W_right.velocity + R.W_right.c * ((gamma+1)/(2*gamma) * (R.P_star/R.W_right.pressure) + (gamma-1)/(2*gamma))**.5

                if(sampling_point > Speed_Shock):
                    return R.W_right
                else:
                    return R.W_right_star

