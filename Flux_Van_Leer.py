class Flux_Van_Leer:

    def __init__(self, W1, W2):

        #Flux +

        Mach = W1.velocity/W1.c

        if (Mach <= -1):
            Mach_flux_1 = 0
            Mach_flux_2 = 0
            Mach_flux_3 = 0
        elif ( -1 < Mach < 1):
            Mach_flux_1 = ((Mach+1)/2)**2
            Mach_flux_2 = ((Mach+1)/2)**2 * ((W1.gamma-1)*Mach+2)
            Mach_flux_3 = ((W1.gamma-1)*W1.velocity + 2 * W1.c)**2 * W1.rho * W1.c * Mach_flux_1 / (2 * (W1.gamma+1) * (W1.gamma-1))
        else:
            Mach_flux_1 = Mach
            Mach_flux_2 = W1.gamma * Mach**2 +1
            Mach_flux_3 = W1.rho * W1.c**3 * Mach * (Mach**2/2 + 1/ (W1.gamma -1))

        f1_plus = W1.rho * W1.c * Mach_flux_1
        f2_plus = W1.rho * W1.c**2/W1.gamma * Mach_flux_2
        f3_plus = Mach_flux_3 


        #Flux -
        
        Mach = W2.velocity/W2.c

        if (Mach <= -1):
            Mach_flux_1 = Mach 
            Mach_flux_2 = W2.gamma * Mach**2 + 1
            Mach_flux_3 = W2.rho * W2.c**3 * Mach * (Mach**2/2 + 1/ (W2.gamma -1))

        elif ( -1 < Mach < 1):
            Mach_flux_1 = -((Mach-1)/2)**2
            Mach_flux_2 = -((Mach-1)/2)**2 * ((W2.gamma-1)*Mach-2)
            Mach_flux_3 = ((W2.gamma-1)*W2.velocity - 2 * W2.c)**2 * W2.rho * W2.c * Mach_flux_1 / (2 * (W2.gamma+1) * (W2.gamma-1))
        else:
            Mach_flux_1 = 0
            Mach_flux_2 = 0
            Mach_flux_3 = 0

        f1_minus = W2.rho * W2.c * Mach_flux_1
        f2_minus = W2.rho * W2.c**2/W2.gamma * Mach_flux_2
        f3_minus = Mach_flux_3 

        self.f1 = f1_plus + f1_minus
        self.f2 = f2_plus + f2_minus
        self.f3 = f3_plus + f3_minus

