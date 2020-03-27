from Flux import Flux
from State import State

#Richtmeyer computation of the predictor
def Predictor(W_L, W_R, delta_t, delta_x):
    F_L = Flux(W_L)
    F_R = Flux(W_R)

    u1 = .5 * (W_R.u1 + W_L.u1) - delta_t/delta_x/2 * (F_R.f1 - F_L.f1)
    u2 = .5 * (W_R.u2 + W_L.u2) - delta_t/delta_x/2 * (F_R.f2 - F_L.f2)
    u3 = .5 * (W_R.u3 + W_L.u3) - delta_t/delta_x/2 * (F_R.f3 - F_L.f3)

    density = u1
    velocity = u2/u1
    pressure = (u3 - u2**2/(2*u1)) * (W_L.gamma - 1)

    return State(None, W_L.gamma, density, velocity, pressure)

