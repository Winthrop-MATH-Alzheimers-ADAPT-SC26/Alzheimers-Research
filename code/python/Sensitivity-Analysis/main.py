from model_computation import ODEModel
from non_dimensionalized import NonDimODEModel
import numpy as np
import pandas as pd

RESULTS_CSV = "sobol_results.csv"

if __name__ == '__main__':
    # given arbitrary parameters
    params = [65641/10000, 15778.463/10000, 315569260, 6311385.2, 52.2958, 1.78467,
            15778.463/10000, 0.07176, 0.398406, 146.308032/10000, 86.84/10000, 199.16/10000,  # c3 originally 0.1, testing with 1.0 for analysis.
            3035.98/10000, 3155692.6, 10.9999/6, 0.3588, 0.8684/10, 100, 100]

    initial_conditions = [0.0, 100, 0.0, 0.0, 0.0] # initial conditions [ABeta_0, Ca_0, Tau_0, N_0, C_0]
    t_span = (0, 50)
    t_eval = np.linspace(*t_span, 100)
    model_no_treatment = ODEModel(params, initial_conditions, t_span, t_eval)

    # start solver for ODEs
    # print("...starting solver...")
    # model_no_treatment.results()

    sample_sizes = 10
    print(f"...starting sensitivity analysis for n={sample_sizes}...")
    sobol_analysis = model_no_treatment.sensitivity_analysis(sample_sizes)
    model_no_treatment.sobol_dataframe_output(sobol_analysis)
