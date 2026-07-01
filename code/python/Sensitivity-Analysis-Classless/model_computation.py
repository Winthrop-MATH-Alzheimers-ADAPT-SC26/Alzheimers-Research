import numpy as np
from scipy.integrate import solve_ivp, odeint
from numba import njit

@njit
def alzheimers_ode_numba(y, t, params):
    Ab, Ca, Tau, N, C = y
    
    a1, a2, b1, b2, c1, c2, c3, d1, d2, e1, e2, e3, k1, k2, k3, k4, k5, sig1, sig2 = params
    
    dAb = a1 + a2 * (Ca / (Ca + sig1)) - k1 * Ab
    dCa = b1 + b2 * Ab - k2 * Ca
    dTau = c1 + c2 * Ab + c3 * (Ca / (Ca + sig2)) - k3 * Tau
    dN = d1 + d2 * Tau - k4 * N
    dC = e1 + e2 * N + e3 * Tau - k5 * C
    
    return (dAb, dCa, dTau, dN, dC)

def solve_for_params(args):
    current_params, init_cond, t_span, t_eval = args
    
    # must be a np array for the solver
    params_array = np.array(current_params, dtype = np.float64)
    
    # run ODE solver (uses LSODA)
    sol = odeint(
        func = alzheimers_ode_numba, 
        y0 = init_cond, 
        t = t_eval, 
        args = (params_array,) 
    )
    
    return sol[-1, 4]