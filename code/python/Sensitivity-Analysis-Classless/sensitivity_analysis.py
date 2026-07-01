from model_computation import run_batch, alzheimers_ode_numba
import numpy as np
import pandas as pd
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol
import os
import time

if __name__ == '__main__':
    params = (65641/10000, 15778.463/10000, 315569260, 6311385.2, 52.2958, 1.78467, 15778.463/10000, 0.07176, 0.398406, 146.308032/10000, 86.84/10000, 199.16/10000, 3035.98/10000, 3155692.6, 10.9999/6, 0.3588, 0.8684/10, 100, 100)

    # get plus/minus 15% bounds for parameters
    lower_bound = [param * 0.85 for param in params]
    upper_bound = [param * 1.15 for param in params]
    total_bounds = list(zip(lower_bound, upper_bound))

    paramNames = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'c3', 'd1', 'd2', 'e1', 'e2', 'e3', 'k1', 'k2', 'k3', 'k4', 'k5', 'sig1', 'sig2']

    problem = {
        'num_vars': len(paramNames),
        'names': paramNames,
        'bounds': total_bounds
    }

    # log_2 (n) sample size
    sample_size = 2**20

    # create args batch and convert to jnp array
    print("...sampling parameter values...")
    param_values = sobol_sample.sample(problem, sample_size, calc_second_order = False)
    param_values = np.ascontiguousarray(param_values, dtype = np.float64)

    u0 = np.array([0.0, 100.0, 0.0, 0.0, 0.0], dtype = np.float64)
    t_span = (0, 50)
    t_eval = np.linspace(t_span[0], t_span[1], 100, dtype = np.float64)

    # multiprocessing ODE solving
    print("...starting CPU solver...")
    num_cores = os.cpu_count()
    print(f"...evaluating {len(param_values)} parameter sets across {num_cores} cores...")

    start_time = time.perf_counter()
    model_out = run_batch(alzheimers_ode_numba.address, param_values, u0, t_eval)
    end_time = time.perf_counter()

    total_time = end_time - start_time
    sims_per_sec = len(param_values) / total_time
    print(f"...batch took {total_time:.2f} seconds with {sims_per_sec:.1f} simulations per second...")

    # analyze results from simulations
    print("...analyzing...")
    sobol_analysis = sobol.analyze(
        problem, 
        model_out, 
        calc_second_order = False,
        n_processors = num_cores,
        print_to_console = False)
    
    # analysis output
    total = sobol_analysis['ST']
    conf = sobol_analysis['ST_conf']

    df = pd.DataFrame({
        'Params': paramNames,
        'Sobol First Order (S1)': sobol_analysis['S1'],
        'Sobol Total (ST)': total,
        'Sobol Total Conf': conf,
        'High Accuracy': np.where((total / conf) < 10, '--', 'Yes'),
        'Below Relative Cutoff': np.where(total < (max(total) * 0.01), 'Yes', 'No')
    }).sort_values(by = 'Sobol Total (ST)', ascending = False)

    print("\n ...Sobol Analysis Table...")
    print(df.to_string(index=False))