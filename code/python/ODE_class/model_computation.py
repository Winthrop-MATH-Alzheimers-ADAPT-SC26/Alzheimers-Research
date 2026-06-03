import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sensitivity_analysis import SensitivityAnalysis
# define the class for system of equations
class ODEModel:
    def __init__(self, params, initial_conditions, t_span, t_eval): 
        # params list is ordered: # a1, a2, b1, b2, c1, c2, c3,
        # d1, d2, e1, e2, e3, u1, u2,
        # e3, u1, u2, u3, k1, k2, k3,
        # k4, k5, sigma, r.
        self.a1 = params[0]
        self.a2 = params[1]
        self.b1 = params[2]
        self.b2 = params[3]
        self.c1 = params[4]
        self.c2 = params[5]
        self.c3 = params[6]
        self.d1 = params[7]
        self.d2 = params[8]
        self.e1 = params[9]
        self.e2 = params[10]
        self.e3 = params[11]
        #self.u1 = params[12] treatment parameters not needed here,
        #self.u2 = params[13] implemented in optimal control later
        #self.u3 = params[14]
        self.k1 = params[12]
        self.k2 = params[13]
        self.k3 = params[14]
        self.k4 = params[15]
        self.k5 = params[16]
        self.sigma = params[17]
        self.r = params[18]
        # for modeling solution space
        self.init_cond = initial_conditions
        self.t_span = t_span
        self.t_eval = t_eval
    def return_string_list(self):
        params_string = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2',
                         'c3', 'd1', 'd2', 'e1', 'e2', 'e3',
                         'k1', 'k2', 'k3', 'k4', 'k5', 'sig',
                         'r']
        num_params = len(params_string)
        return [params_string, num_params]
    # use call to allow for solutions to be computed later
    def __call__(self, t, y):
        ABeta, Ca, Tau, N, C = y
        dABeta = self.a1 + self.a2 * (Ca**2 / (Ca**2 + self.sigma**2)) - self.k1 * ABeta #- self.u1 * ABeta
        dCa = self.b1 + self.b2 * ABeta - self.k2 * Ca # - self.u2 * Ca
        dTau = self.c1 + self.c2 * ABeta + self.c3 * Ca - self.k3 * Tau #- self.u3 * Tau
        dN = self.d1 + self.d2 * Tau - self.k4 * N
        dC = self.e1 + self.e2 * N * self.r + self.e3 * Tau - self.k5 * C
        return [dABeta, dCa, dTau, dN, dC]
    # call sensitivity analysis using model in class to run SAlib sobol analysis
    def sensitivity_analysis(self):
        print(SensitivityAnalysis(self).run_sensitivity_analysis())
    def solution(self):
        solution = solve_ivp(self, self.t_span, self.init_cond, t_eval=self.t_eval, method='LSODA')
        return solution
    def results(self):
        print("...calculated iterations...")
        print(*self.solution().y, sep=", ")

# given arbitrary parameters
params = [65641, 15778.463, 315569260, 6311385.2, 52.2958, 1.78467,
          3.62493, 0.07176, 0.0195, 146.308032, 86.84, 199.16, 0.5, 0.5, 0.5,
          3035.98, 21874000, 10.9999, 0.003588, 0.8684, 0.00027, 1]
initial_conditions = [0.0, 0.0, 0.0, 0.0, 0.0] # initial conditions [ABeta_0, Ca_0], 0.0, 0.0, 0.0
t_span = (0, 100)
t_eval = np.linspace(*t_span, 100)
model_no_treatment = ODEModel(params, initial_conditions, t_span, t_eval)
# start solver for ODEs
print("...starting solver...")
model_no_treatment.results()

# create a grid with 2 rows and 2 columns
#fig, axs = plt.subplots(2, 2)
# access a specific subplot using 2D indexing [row, col]
#axs[0, 0].plot(t_eval, sol.y[0])
#axs[0, 1].plot(t_eval, sol.y[1])
#plt.show()

