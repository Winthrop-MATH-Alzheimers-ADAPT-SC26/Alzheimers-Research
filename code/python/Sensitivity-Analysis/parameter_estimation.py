import numpy as np
import pandas as pd

class ParameterEstimation:
    def __init__(self, params): # have path as arg when we get dataset
        #self.df = pd.read_csv(path)
        self.params = params
    
    # first test with -15%, +15% of our estimated params from the literature
    # using df, get estimated param ranges
    def get_param_ranges(self):
        lower_bound = [param * 0.85 for param in self.params] # multiply 0.85 to subtract 15%
        upper_bound = [param * 1.15 for param in self.params]
        total_bounds = list(zip(lower_bound, upper_bound))
        return total_bounds
