import numpy as np
import pandas as pd

class ParameterEstimation:
    def __init__(self, path):
        self.df = pd.read_csv(path)

    # using df, get estimated param ranges
    def get_param_ranges(self):
        # calculate range for every column in the dataframe at once
        all_ranges = self.df.apply(lambda x: x.max() - x.min())
        return all_ranges

