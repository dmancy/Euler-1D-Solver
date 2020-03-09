from Flux import Flux
import numpy as np

class Flux_Roe:

    def __init__(self, W_L, W_R):

        F_L = Flux(W_L)
        F_R = Flux(W_R)

        #Roe variables
        u_tilde = (np.sqrt(W_L.rho) * W_L.velocity + np.sqrt(W_R.rho) * W_R.velocity)/(np.sqrt(W_L.rho) + np.sqrt(W_R.rho))

        h_L = W_L.pressure /((W_L.gamma-1) * W_L.rho) + W_L.pressure / W_L.rho 
        h_R = W_R.pressure /((W_R.gamma-1) * W_R.rho) + W_R.pressure / W_R.rho 

        h_tilde = (np.sqrt(W_L.rho) * h_L + np.sqrt(W_R.rho) * h_R)/(np.sqrt(W_L.rho) + np.sqrt(W_R.rho))


        c_tilde = np.sqrt((W_L.gamma - 1) * (h_tilde - .5 * u_tilde**2))

        M_tilde = u_tilde/c_tilde


        #Roe Eigenvalues

        lambda_1 = abs(u_tilde - c_tilde)
        lambda_2 = abs(u_tilde)
        lambda_3 = abs(u_tilde + c_tilde)


        delta_rho = (W_R.u1 - W_L.u1)
        delta_momentum = (W_R.u2 - W_L.u2)
        delta_rho_E = (W_R.u3 - W_L.u3)

        alpha_1 = 0.25 * M_tilde * (2 + (W_L.gamma - 1) * M_tilde) * delta_rho \
                - .5/c_tilde * (1 + (W_L.gamma - 1) * M_tilde) * delta_momentum \
                + (W_L.gamma - 1)/2 * 1/c_tilde**2 * delta_rho_E 

        alpha_2 = (1 - (W_L.gamma - 1)/2 * M_tilde**2) * delta_rho \
                + (W_L.gamma -1) * M_tilde/c_tilde * delta_momentum\
                - (W_L.gamma - 1) * delta_rho_E/c_tilde**2
        
        alpha_3 = -.25 * M_tilde * (2 - (W_L.gamma - 1) * M_tilde) * delta_rho \
                + 1/(2*c_tilde) * (1- ( W_L.gamma - 1) * M_tilde) * delta_momentum \
                + (W_L.gamma - 1)/2 * delta_rho_E/c_tilde**2

        self.f1 = .5 * (F_R.f1 + F_L.f1) -.5 * (lambda_1 * alpha_1 + lambda_2 * alpha_2 + lambda_3 * alpha_3)

        self.f2 = .5 * (F_R.f2 + F_L.f2) - .5 * (lambda_1 * alpha_1 * (u_tilde - c_tilde) \
                + lambda_2 * alpha_2 * u_tilde + lambda_3 * alpha_3 * (u_tilde + c_tilde))

        self.f3 = 0.5 * (F_R.f3 + F_L.f3) - .5 * (lambda_1 * alpha_1 * (h_tilde - u_tilde * c_tilde) \
                + .5 * lambda_2 * alpha_2 * u_tilde**2 + lambda_3 * alpha_3 * (h_tilde + u_tilde * c_tilde))


