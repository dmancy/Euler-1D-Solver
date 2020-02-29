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

    def print(self):
        """Print density, velocity, pressure, temperature and sound speed in this region"""
        print("State : {} \nDensity : {} \nVelocity : {} \nPressure : {} \nTemperature : {} \nSound speed : {}\n".format(self.name, self.rho, self.velocity, self.pressure, self.T, self.c))
