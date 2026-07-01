import numpy as np
from numba import cfunc, carray
from numbalsoda import lsoda_sig, lsoda

@cfunc(lsoda_sig)
def alzheimers_ode_numba(t, u, du, p):
    # cast pointers into arrays
    y = carray(u, (5,))
    dydt = carray(du, (5,))
    params = carray(p, (19,))

    Ab, Ca, Tau, N, C = y[0], y[1], y[2], y[3], y[4]
    
    dydt[0] = params[0] + params[1] * (Ca / (Ca + params[17])) - params[12] * Ab
    dydt[1] = params[2] + params[3] * Ab - params[13] * Ca
    dydt[2] = params[4] + params[5] * Ab + params[6] * (Ca / (Ca + params[18])) - params[14] * Tau
    dydt[3] = params[7] + params[8] * Tau - params[15] * N
    dydt[4] = params[9] + params[10] * N + params[11] * Tau - params[16] * C

def solve_for_params(args):
    current_params, init_cond, t_span, t_eval = args
    
    # must be a np array for the solver
    params_array = np.array(current_params, dtype = np.float64)
    u0 = np.array(init_cond, dtype = np.float64)
    
    usol, success = lsoda(
        funcptr = alzheimers_ode_numba.address,
        u0 = u0,
        t_eval = t_eval,
        data = params_array
    )
    
    return usol[-1, 4]