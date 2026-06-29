from model_computation import ODEModel
from non_dimensionalized import NonDimODEModel
import numpy as np
import pandas as pd

RESULTS_CSV = "sobol_results.csv"

if __name__ == '__main__':
    # given arbitrary parameters
    params = [65641/10000, 15778.463/10000, 315569260, 6311385.2, 52.2958, 1.78467,
            15778.463/10000, 0.07176, 0.398406, 146.308032/10000, 86.84/10000, 199.16/10000,  # c3 originally 0.1, testing with 1.0 for analysis.
            3035.98/10000, 3155692.6, 10.9999/6, 0.3588, 0.8684/10, 100, 100, 1]
    def reparam(params):
        a1 = params[0]
        a2 = params[1]
        b1 = params[2]
        b2 = params[3]
        c1 = params[4]
        c2 = params[5]
        c3 = params[6]
        d1 = params[7]
        d2 = params[8]
        e1 = params[9]
        e2 = params[10]
        e3 = params[11]
        k1 = params[12]
        k2 = params[13]
        k3 = params[14]
        k4 = params[15]
        k5 = params[16]
        sigma = params[17]
        sigma2 = params[18]
        r = params[19]

        a1_bar = a2 / a1
        sigma1_bar = sigma * k1 / b1
        k2_bar = k2 / k1
        b2_bar = b2 * a1 / (b1 * k1)
        k3_bar = k3 / k1
        c2_bar = c2 * a1 / (c1 * k1)
        c3_bar = c3 / c1
        sigma2_bar = sigma2 * k1 / b1
        k4_bar = k4 / k1
        d2_bar = d2 * c1 / (d1 * k1)
        k5_bar = k5 / k1
        e2_bar = e2 * d1 / (e1 * k1)
        e3_bar = e3 * c1 / (e1 * k1)

        return [a1_bar, sigma1_bar, k2_bar, b2_bar, k3_bar, c2_bar, c3_bar,
                sigma2_bar, k4_bar, d2_bar, k5_bar, e2_bar, e3_bar]

    initial_conditions = [0.0, (100 * (3035.98/10000))/315569260, 0.0, 0.0, 0.0] # initial conditions [ABeta_0, Ca_0, Tau_0, N_0, C_0]
    t_span = (0, 50 * (3035.98/10000))
    t_eval = np.linspace(*t_span, 100)
    model_no_treatment = NonDimODEModel(reparam(params), initial_conditions, t_span, t_eval)

    # start solver for ODEs
    # print("...starting solver...")
    # model_no_treatment.results()

    sample_sizes = 18
    print(f"...starting sensitivity analysis for n={sample_sizes}...")
    sobol_analysis = model_no_treatment.sensitivity_analysis(sample_sizes)
    model_no_treatment.sobol_dataframe_output(sobol_analysis)
