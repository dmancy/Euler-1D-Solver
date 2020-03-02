class Flux_Zha_Bilgen:

    def __init__(self, W1, W2):

        #Flux +
        
        U = max(W1.velocity, 0)
        Mach = W1.velocity/W1.c

        if (Mach <= -1):
            p_plus = 0
            pu_plus = 0
        elif (-1 < Mach < 1):
            p_plus = W1.pressure * (Mach +1)/2
            pu_plus = W1.pressure * (W1.velocity + W1.c)/2
        else:
            p_plus = W1.pressure
            pu_plus = W1.pressure * W1.velocity

        f1_plus = U * W1.rho
        f2_plus = U * W1.rho * W1.velocity + p_plus
        f3_plus = U * (W1.pressure/(W1.gamma-1) + W1.rho * W1.velocity**2/2) + pu_plus
        

        #Flux -
        
        U = min(W2.velocity, 0)
        Mach = W2.velocity/W2.c

        if (Mach <= -1):
            p_minus = W2.pressure
            pu_minus = W2.pressure * W2.velocity
        elif (-1 < Mach < 1):
            p_minus = W2.pressure * (1 - Mach)/2
            pu_minus = W2.pressure * (W2.velocity - W2.c)/2
        else:
            p_minus = 0
            pu_minus = 0

        f1_minus = U * W2.rho
        f2_minus = U * W2.rho * W2.velocity + p_minus
        f3_minus = U * (W2.pressure/(W2.gamma-1) + W2.rho * W2.velocity**2/2) + pu_minus


        #Total Flux
        self.f1 = f1_plus + f1_minus
        self.f2 = f2_plus + f2_minus
        self.f3 = f3_plus + f3_minus
