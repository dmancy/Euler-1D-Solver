import math

class State:
    """Define the state of an area by:
    - name
    - gamma
    - density rho
    - velocity u
    - pressure p
    - Temperature T
    - sound speed c"""

    def __init__(self, title, gam, density, velocity, pressure):
        """Constructor"""
        self.name = title
        self.gamma = gam
        self.R = 287.058
        self.rho = density
        self.velocity = velocity
        self.pressure = pressure
        self.T = self.pressure/(self.R * self.rho)
        self.c = math.sqrt(self.gamma * self.R * self.T)

        #Conservative variables
        self.u1 = self.rho
        self.u2 = self.rho * self.velocity
        self.u3 = self.pressure/(self.gamma - 1) + self.rho * self.velocity**2/2

    def print(self):
        """Print density, velocity, pressure, temperature and sound speed in this region"""
        print("State : {} \nDensity : {} \nVelocity : {} \nPressure : {} \nTemperature : {} \nSound speed : {}\n".format(self.name, self.rho, self.velocity, self.pressure, self.T, self.c))

    def update_primitive_variables(self):
        """ Update primitive variables from the current conservative variables"""

        self.rho = self.u1
        self.velocity = self.u2/self.u1
        self.pressure = (self.u3 - self.u2**2/(2*self.u1)) * (self.gamma - 1)
        self.T = self.pressure/(self.R * self.rho)
        self.c = math.sqrt(self.gamma * self.R * self.T)

    def solver_step(self, delta_t, delta_x, F_L, F_R):
        """ Update the state vector with an explicit scheme Ui = Ui - delta_t/delta_x * (F_R - F_L) """
        self.u1 -= delta_t / delta_x * (F_R.f1 - F_L.f1)
        self.u2 -= delta_t / delta_x * (F_R.f2 - F_L.f2)
        self.u3 -= delta_t / delta_x * (F_R.f3 - F_L.f3)

        self.update_primitive_variables()


