import numpy as np
from numba import cfunc, carray, njit, prange, types
from numbalsoda import lsoda_sig, lsoda



@cfunc(lsoda_sig)
def alzheimers_ode_numba(t, u, du, p):
    # cast pointers into arrays
    y = carray(u, (5,))
    dydt = carray(du, (5,))
    params = carray(p, (19,))

    Ab, Ca, Tau, N, C = y[0], y[1], y[2], y[3], y[4]
    
    a1, a2, b1, b2, c1, c2, c3, d1, d2, e1, e2, e3, k1, k2, k3, k4, k5, sig1, sig2 = params

    dydt[0] = a1 + a2 * (Ca / (Ca + sig1)) - k1 * Ab 
    dydt[1] = b1 + b2 * Ab - k2 * Ca 
    dydt[2] = c1 + c2 * Ab + c3 * (Ca / (Ca + sig2)) - k3 * Tau 
    dydt[3] = d1 + d2 * Tau - k4 * N
    dydt[4] = e1 + e2 * N + e3 * Tau - k5 * C



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