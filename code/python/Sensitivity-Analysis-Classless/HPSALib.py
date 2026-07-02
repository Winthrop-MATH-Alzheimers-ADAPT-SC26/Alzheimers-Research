# This code is adapted from sobol.py from SALib.
from numba import njit, prange
import numpy as np
from scipy.stats import qmc, norm



# Numba optimized Saltelli sampling (first order only)
@njit(parallel = True)
def _build_saltelli_first_order(base_sequence):
    N = base_sequence.shape[0]
    D = base_sequence.shape[1] // 2

    # initialize empty array
    saltelli_sequence = np.empty((N * (D + 2), D), dtype = base_sequence.dtype)

    # this loop is parallelized by numba
    for i in prange(N):
        start_idx = i * (D + 2)

        for j in range(D):
            saltelli_sequence[start_idx, j] = base_sequence[i, j]

        for k in range(D):
            row_idx = start_idx + 1 + k

            for j in range(D):
                if j == k:
                    saltelli_sequence[row_idx, j] = base_sequence[i, j + D]
                else:
                    saltelli_sequence[row_idx, j] = base_sequence[i, j]

        last_row_idx = start_idx + 1 + D

        for j in range(D):
            saltelli_sequence[last_row_idx, j] = base_sequence[i, j + D]

    return saltelli_sequence

def saltelli_sample(problem, N):
    # generate the base sequence
    qrng = qmc.Sobol(d = 2 * problem["num_vars"])
    base_sequence = qrng.random(N)
    param_values = _build_saltelli_first_order(base_sequence)

    # scale the sobol sequence to our parameter bounds
    bounds = np.array(problem["bounds"])
    lower_bounds = bounds[:, 0]
    ranges = bounds[:, 1] - bounds[:, 0]
    param_values = param_values * ranges + lower_bounds

    return param_values



# Numba optimized Sobol analysis (first order only)
@njit(parallel = True)
def _numba_sobol_kernel(A, B, AB, r, Z):
    N, D = AB.shape
    num_resamples = r.shape[0]
    
    S1 = np.zeros(D)
    ST = np.zeros(D)
    S1_conf = np.zeros(D)
    ST_conf = np.zeros(D)
    
    # precompute the variance of the base array
    sum_AB = 0.0
    for i in range(N):
        sum_AB += A[i] + B[i]
    mean_y = sum_AB / (2.0 * N)
    
    var_y_num = 0.0
    for i in range(N):
        var_y_num += (A[i] - mean_y)**2 + (B[i] - mean_y)**2
    var_y = var_y_num / (2.0 * N)
    
    # numba parallel loop
    for j in prange(D):
        # base sensitivity computations
        s1_num = 0.0
        st_num = 0.0
        for i in range(N):
            s1_num += B[i] * (AB[i, j] - A[i])
            st_num += (A[i] - AB[i, j])**2
            
        S1[j] = (s1_num / N) / var_y
        ST[j] = (0.5 * st_num / N) / var_y
        
        # inline bootstrap resampling
        s1_resamples = np.empty(num_resamples)
        st_resamples = np.empty(num_resamples)
        
        for b in range(num_resamples):
            sum_y_r = 0.0
            # compute resampled mean
            for i in range(N):
                idx = r[b, i]
                sum_y_r += A[idx] + B[idx]
            mean_y_r = sum_y_r / (2.0 * N)
            
            var_y_r_num = 0.0
            s1_num_r = 0.0
            st_num_r = 0.0
            
            # compute resampled variances and estimators
            for i in range(N):
                idx = r[b, i]
                a_val = A[idx]
                b_val = B[idx]
                ab_val = AB[idx, j]
                
                var_y_r_num += (a_val - mean_y_r)**2 + (b_val - mean_y_r)**2
                s1_num_r += b_val * (ab_val - a_val)
                st_num_r += (a_val - ab_val)**2
                
            var_y_r = var_y_r_num / (2.0 * N)
            
            s1_resamples[b] = (s1_num_r / N) / var_y_r
            st_resamples[b] = (0.5 * st_num_r / N) / var_y_r
            
        # calculate sample SD of resamples
        s1_r_mean = 0.0
        st_r_mean = 0.0
        for b in range(num_resamples):
            s1_r_mean += s1_resamples[b]
            st_r_mean += st_resamples[b]
        s1_r_mean /= num_resamples
        st_r_mean /= num_resamples
        
        s1_r_var = 0.0
        st_r_var = 0.0
        for b in range(num_resamples):
            s1_r_var += (s1_resamples[b] - s1_r_mean)**2
            st_r_var += (st_resamples[b] - st_r_mean)**2
            
        S1_conf[j] = Z * np.sqrt(s1_r_var / (num_resamples - 1))
        ST_conf[j] = Z * np.sqrt(st_r_var / (num_resamples - 1))
        
    return S1, S1_conf, ST, ST_conf

def fast_sobol_analyze(problem, Y, num_resamples = 100, conf_level = 0.95, seed = None):
    D = problem["num_vars"]
    N = int(Y.size / (D + 2))
    
    # normalize the model output matching SALib's scaling step
    Y_norm = (Y - Y.mean()) / Y.std()
    
    # fast unrolled slicing instead of a slow loop
    step = D + 2
    A = Y_norm[0 : Y_norm.size : step]
    B = Y_norm[step - 1 : Y_norm.size : step]
    
    AB = np.empty((N, D), dtype = Y.dtype)
    for j in range(D):
        AB[:, j] = Y_norm[(j + 1) : Y_norm.size : step]
        
    # generate the bootstrap index grid
    if seed:
        rng = np.random.default_rng(seed)
    else:
        rng = np.random.default_rng()
        
    # shape (num_resamples, N) allows Numba to read memory sequentially along rows inside the inner loops
    r = rng.integers(0, N, size = (num_resamples, N))
    Z = norm.ppf(0.5 + conf_level / 2)
    
    S1, S1_conf, ST, ST_conf = _numba_sobol_kernel(A, B, AB, r, Z)
    
    return {
        "S1": S1,
        "S1_conf": S1_conf,
        "ST": ST,
        "ST_conf": ST_conf
    }