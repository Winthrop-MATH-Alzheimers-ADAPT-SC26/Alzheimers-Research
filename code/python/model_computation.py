import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class ODEModel:
    def __init__(self, a1, a2, b1, b2, c1, c2, c3, d1, e1,e2, k1, k2, k3, alph1, alph2, alph3, delt1, delt2, delt3, sigma, kn, kc):
        self.a1 = a1
        self.a2 = a2
        self.b1 = b1
        self.b2 = b2
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.d1 = d1
        self.e1 = e1
        self.e2 = e2
        self.sigma = sigma
        self.kn = kn
        self.kc = kc
        self.r1 = k1 + (1 - alph1) * delt1
        self.r2 = k2 + (1 - alph2) * delt2
        self.r3 = k3 + (1 - alph3) * delt3

    def __call__(self, t, y):
        ABeta, Ca, Tau, N, C = y
        dABeta = self.a1 - self.r1 * ABeta + self.a2 * Ca**2 / (Ca**2 + self.sigma**2)
        dCa = self.b1 - self.r2 * Ca + self.b2 * ABeta
        dTau = self.c1 - self.r3 * Tau + self.c2 * ABeta + self.c3 * Ca
        dN = self.d1 * Tau * (1 - (N/self.kn))
        dC = (self.e1 * N + self.e2 * Tau) * (1 - (C / self.kc))
        return [dABeta, dCa, dTau, dN, dC]



# given arbitrary parameters
#model_eff_50 = ODEModel(a1=0.01, a2=0.05, b1=1000, b2=20,
#                 k1=0.00002617, k2=0.1, alph1=0.5, alph2=0.5,
#                 delt1=0.4, delt2=0.4, sigma=0.00027)

y0 = [0.0, 0.0] # initial conditions [ABeta_0, Ca_0]
t_span = (0, 100)
t_eval = np.linspace(*t_span, 1000)

#sol = solve_ivp(model_eff_50, t_span, y0, t_eval=t_eval, method='RK45')

#sol.y[0] #→ ABeta over time
#sol.y[1] #→ Ca over time

#print(sol.y[0])
#print(sol.y[1])

# create a grid with 2 rows and 2 columns
#fig, axs = plt.subplots(2, 2)

# access a specific subplot using 2D indexing [row, col]
#axs[0, 0].plot(t_eval, sol.y[0])
#axs[0, 1].plot(t_eval, sol.y[1])
#plt.show()

