import numpy as np
from numba import cfunc, carray, njit, prange
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

@njit(parallel = True)
def run_batch(funcptr, param_matrix, u0, t_eval):
    n_sims = param_matrix.shape[0]
    results = np.empty(n_sims, dtype = np.float64)

    for i in prange(n_sims):
        usol, success = lsoda(
            funcptr = funcptr,
            u0 = u0,
            t_eval = t_eval,
            data = param_matrix[i]
        )

        results[i] = usol[-1, 4]

    return results