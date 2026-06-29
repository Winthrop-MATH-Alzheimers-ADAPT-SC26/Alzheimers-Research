import os
import pandas as pd
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol
from vectorized_solver import solve_batch
from tqdm import tqdm

params = [6.5641, 1.5778, 315569260, 6311385.2, 52.2958, 1.78367, 1.5778, 0.0716, 0.398406, 0.01463, 0.008684, 0.01992, 0.303598, 3155692.6, 1.8333, 0.3588, 0.08684, 100, 100]
y0 = jnp.array([0, 100, 0, 0, 0], dtype = jnp.float64)

# get plus/minus 15% bounds for parameters
lower_bound = [param * 0.85 for param in params]
upper_bound = [param * 1.15 for param in params]
total_bounds = list(zip(lower_bound, upper_bound))

paramNames = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'c3', 'd1', 'd2', 'e1', 'e2', 'e3', 'k1', 'k2', 'k3', 'k4', 'k5', 'sig', 'sig2']

problem = {
    'num_vars': len(paramNames),
    'names': paramNames,
    'bounds': total_bounds
}

# log_2 (n) sample size
sample_size = 2**14

# create args batch and convert to jnp array
print("...sampling parameter values...")
param_values = sobol_sample.sample(problem, sample_size, calc_second_order = False)
jnp_param_values = jnp.asarray(param_values)

# run gpu solver
print("...starting gpu solver...")
chunk_size = 16384
num_samples = len(jnp_param_values)
all_results = []

# tqdm progress bar
for i in tqdm(range(0, num_samples, chunk_size), desc = "GPU Batch Processing", unit = "chunk"):
    # get param chunk
    param_chunk = jnp_param_values[i : i + chunk_size]

    # solve chunk of params
    solution_object = solve_batch(y0, param_chunk)
    
    # get final c values from solver
    final_C_values = solution_object.ys[:, -1, 4]
    
    # cpu wait until gpu calculations are finished
    final_C_values.block_until_ready()
    
    # add results to array
    all_results.append(np.array(final_C_values))

# turn into 1d array for results
model_out = np.concatenate(all_results)

# analyze results
print("...analyzing results...")
sobol_analysis = sobol.analyze(
    problem,
    model_out,
    calc_second_order = False,
    n_processors = os.cpu_count()
)

# print to console
total = sobol_analysis['ST']
conf = sobol_analysis['ST_conf']

df = pd.DataFrame({
    'Params': paramNames,
    'Sobol First Order (S1)': sobol_analysis['S1'],
    'Sobol Total (ST)': total,
    'Sobol Total Conf': conf,
    'High Accuracy': np.where((total / conf) < 10, '--', 'Yes'),
    'Below Relative Cutoff': np.where(total < (max(total)*0.01), 'Yes', 'No')
    }).sort_values(by='Sobol Total (ST)', ascending=False)
print("\n ...Sobol Analysis Table...")
print(df.to_string(index=False))
