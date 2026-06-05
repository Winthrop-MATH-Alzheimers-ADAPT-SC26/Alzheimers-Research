from model_computation import ODEModel
import numpy as np
# given arbitrary parameters
params = [65641/10000, 15778.463/10000, 315569260, 6311385.2, 52.2958, 1.78467,
          0.1, 0.07176, 0.398406, 146.308032/10000, 86.84/10000, 199.16/10000,
          3035.98/10000, 3155692.6, 10.9999/6, 0.3588, 0.8684/10, 0.00027, 1]
initial_conditions = [0.0, 200, 0.0, 0.0, 0.0] # initial conditions [ABeta_0, Ca_0, Tau_0, N_0, C_0]
t_span = (0, 100)
t_eval = np.linspace(*t_span, 100)
model_no_treatment = ODEModel(params, initial_conditions, t_span, t_eval)
# start solver for ODEs
print("...starting solver...")
model_no_treatment.results()

print("...starting sensitivity analysis...")
sobol_analysis = model_no_treatment.sensitivity_analysis()
model_no_treatment.sobol_dataframe_output(sobol_analysis)
#model_no_treatment.visualization()
