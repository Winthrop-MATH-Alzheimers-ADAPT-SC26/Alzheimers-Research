import numpy as np
from scipy.integrate import solve_ivp

class ODEModel:
    def __init__(self, a1, a2, b1, b2, k1, k2, alph1, alph2, delt1, delt2, sigma):
        self.a1 = a1
        self.a2 = a2
        self.b1 = b1
        self.b2 = b2
        self.sigma = sigma
        self.r1 = k1 + (1 - alph1) * delt1
        self.r2 = k2 + (1 - alph2) * delt2

    def __call__(self, t, y):
        ABeta, Ca = y
        dABeta = self.a1 - self.r1 * ABeta + self.a2 * Ca**2 / (Ca**2 + self.sigma**2)
        dCa    = self.b1 - self.r2 * Ca + self.b2 * ABeta
        return [dABeta, dCa]


# given arbitrary parameters
model_eff_50 = ODEModel(a1=0.01, a2=0.05, b1=1000, b2=20,
                 k1=0.00002617, k2=0.1, alph1=0.5, alph2=0.5,
                 delt1=0.4, delt2=0.4, sigma=0.00027)

y0 = [0.0, 0.0] # initial conditions [ABeta_0, Ca_0]
t_span = (0, 100)
t_eval = np.linspace(*t_span, 1000)

sol = solve_ivp(model_eff_50, t_span, y0, t_eval=t_eval, method='RK45')

sol.y[0] #→ ABeta over time
sol.y[1] #→ Ca over time

print(sol.y[0])
print(sol.y[1])


