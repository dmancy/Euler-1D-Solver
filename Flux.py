
class Flux:

    def __init__(self, W_state):
        """Flux Constructor"""
        self.f1 = W_state.rho * W_state.velocity
        self.f2 = W_state.rho * W_state.velocity**2 + W_state.pressure
        self.f3 = (W_state.pressure * W_state.gamma / (W_state.gamma-1) + W_state.rho * W_state.velocity**2/2) * W_state.velocity;
