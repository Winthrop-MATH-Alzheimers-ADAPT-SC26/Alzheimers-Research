import numpy as np
import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
from model_computation import ODEModel
class SensitivityAnalysis:
    def __init__(self, model, param_ranges):
        self.model = model
        self.solutions = model.solution().y
        self.param_range = param_ranges
        self.problem = {
            'num_vars': self.model.return_string_list()[1],
            'names': self.model.return_string_list()[0],
            'bounds': self.param_range}

    # using SAlib, perform sensitivity analysis on model
    def run_sensitivity_analysis(self):
        param_values = saltelli.sample(self.problem, 1024)
        # run model 
        Y = self.solutions
        return sobol.analyze(self.problem, Y, print_to_console=True)
