import numpy as np
import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
from parameter_estimation import ParameterEstimation
class SensitivityAnalysis:
    def __init__(self, model):
        self.model = model
        self.param_range = ParameterEstimation('../').get_param_ranges()
        self.problem = {
            'num_vars': self.model.return_string_list()[1],
            'names': self.model.return_string_list()[0],
            'bounds': self.param_range}

    # using SAlib, perform sensitivity analysis on model
    def run_sensitivity_analysis(self):
        param_values = saltelli.sample(self.problem, 1024)
        # run model 
        Y = self.model.evaluate(param_values)
        return sobol.analyze(self.problem, Y, print_to_console=True)
