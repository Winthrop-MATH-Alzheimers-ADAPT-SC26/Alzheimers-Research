from model_computation import ODEModel
import numpy as np
# given arbitrary parameters
params = [65641, 15778.463, 315569260, 6311385.2, 52.2958, 1.78467,
          3.62493, 0.07176, 0.0195, 146.308032, 86.84, 199.16,
          3035.98, 21874000, 10.9999, 0.003588, 0.8684, 0.00027, 1]
initial_conditions = [0.0, 0.0, 0.0, 0.0, 0.0] # initial conditions [ABeta_0, Ca_0, Tau_0, N_0, C_0]
t_span = (0, 100)
t_eval = np.linspace(*t_span, 100)
model_no_treatment = ODEModel(params, initial_conditions, t_span, t_eval)
# start solver for ODEs
print("...starting solver...")
model_no_treatment.results()

model_no_treatment.visualization()
