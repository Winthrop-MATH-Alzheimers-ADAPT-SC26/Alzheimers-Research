import numpy as np
import pandas as pd

class ParameterEstimation:
    def __init__(self, path):
        self.df = pd.read_csv(path)

    # using df, get estimated param ranges
    def get_param_ranges(self):
        all_param_ranges = []
        # for param estimates from literature, hard code param range
        calcium = []

        # calculate range for every column in the dataframe at once
        all_df_ranges = self.df.apply(lambda x: x.max() - x.min())
        all_param_ranges.append(all_df_ranges)
        return all_param_ranges

