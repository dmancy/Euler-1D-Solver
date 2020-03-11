

def Flux_Splitting_Computation(W, lambda1, lambda2, lambda3):

        f1 = (W.gamma - 1)/ W.gamma * W.rho * lambda1 + W.rho/(2*W.gamma) * lambda2 + W.rho/(2*W.gamma) * lambda3
        f2 = (W.gamma - 1)/ W.gamma * W.rho * lambda1 * W.velocity + W.rho/(2*W.gamma) * lambda2 * (W.velocity + W.c) + W.rho/(2*W.gamma) * lambda3 * (W.velocity - W.c)
        f3 = (W.gamma - 1)/ W.gamma * W.rho * lambda1 * W.velocity**2/2 + W.rho/(2*W.gamma) * lambda2 * (W.velocity**2/2 + W.c**2/(W.gamma-1) + W.velocity * W.c) + W.rho/(2*W.gamma) * lambda3 * (W.velocity**2/2 + W.c**2/(W.gamma-1) - W.velocity * W.c)

        return f1, f2, f3




class Flux_Steger_Warming:

    def __init__(self, W_1, W_2):
        """Flux Constructor"""

        #F_plus
        lambda_1 = W_1.velocity
        lambda_1_plus = max(lambda_1, 0)

        lambda_2 = W_1.velocity + W_1.c
        lambda_2_plus = max(lambda_2, 0)

        lambda_3 = W_1.velocity - W_1.c
        lambda_3_plus = max(lambda_3, 0)

        F_plus = Flux_Splitting_Computation(W_1, lambda_1_plus, lambda_2_plus, lambda_3_plus)


        #F_minus
        lambda_1 = W_2.velocity
        lambda_1_minus = min(lambda_1, 0)

        lambda_2 = W_2.velocity + W_2.c
        lambda_2_minus = min(lambda_2, 0)

        lambda_3 = W_2.velocity - W_2.c
        lambda_3_minus = min(lambda_3, 0)

        F_minus = Flux_Splitting_Computation(W_2, lambda_1_minus, lambda_2_minus, lambda_3_minus)

        self.f1 = F_plus[0] + F_minus[0]
        self.f2 = F_plus[1] + F_minus[1]
        self.f3 = F_plus[2] + F_minus[2]

