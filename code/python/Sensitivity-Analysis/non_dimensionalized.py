import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
from parameter_estimation import ParameterEstimation
from SALib.sample import saltelli
from SALib.analyze import sobol
import multiprocessing as mp
from tqdm import tqdm

# outside function for multiprocessing 
def _solve_for_params(args):
    current_params, init_cond, t_span, t_eval = args
    temp_model = NonDimODEModel(current_params, init_cond, t_span, t_eval)
    sol = temp_model.solution()
    return sol.y[4, -1]

# define the class for system of equations
class NonDimODEModel:
    def __init__(self, params, initial_conditions, t_span, t_eval): 
        """initialization.

        Args:
            params: use list of parameter values obtained from literature. Ordered 
            [a1, a2, b1, b2, c1, c2, c3, d1, d2, e1, e2, e3, k1, k2, k3, k4, k5, sigma, r].

            initial_conditions: list of intial condition for state variables.

            t_span: linspace for sequence of time.

            t_eval: use t_span with subdivisions for evaluation.
        ---
        Returns:
            None.
        """
        # for modeling solution space
        self.init_cond = initial_conditions
        self.t_span = t_span
        self.t_eval = t_eval


        self.a1_bar = params[0]
        self.sigma1_bar = params[1]
        self.k2_bar = params[2]
        self.b2_bar = params[3]
        self.k3_bar = params[4]
        self.c2_bar = params[5]
        self.c3_bar = params[6]
        self.sigma2_bar = params[7]
        self.k4_bar = params[8]
        self.d2_bar = params[9]
        self.k5_bar = params[10]
        self.e2_bar = params[11]
        self.e3_bar = params[12]

        self.param_ranges = ParameterEstimation(params).get_param_ranges()

    # used in sobol sensitivity analysis for problem dict
    def return_string_list(self):
        """Using parameter strings, return list of parameter names and length of list."""
        params_string = ['a1_bar', 'sigma1_bar', 'k2_bar', 'b2_bar', 'k3_bar', 'c2_bar', 'c3_bar',
                        'sigma2_bar', 'k4_bar', 'd2_bar', 'k5_bar', 'e2_bar', 'e3_bar'] #, 'u1', 'u2', 'u3'
        num_params = len(params_string)
        return [params_string, num_params]

    # non-dimensionalized system
    def __call__(self, t, y):
        """Use this throughout to call ODE model; useful for using self keyword within class methods."""
        abeta_bar, ca_bar, tau_bar, n_bar, c_bar = y
        dABeta = 1 - abeta_bar - ((self.a1_bar * ca_bar) / (ca_bar + self.sigma1_bar)) - abeta_bar
        dCa = 1 - self.k2_bar * ca_bar + self.b2_bar * abeta_bar - ca_bar
        dTau = 1 - self.k3_bar * tau_bar + self.c2_bar * abeta_bar + ((self.c3_bar * ca_bar) / (ca_bar + self.sigma2_bar)) - tau_bar
        dN = 1 - self.k4_bar * n_bar - self.d2_bar * tau_bar
        dC = 1 - self.k5_bar * c_bar + self.e2_bar * n_bar + self.e3_bar * tau_bar
        return [dABeta, dCa, dTau, dN, dC]

    # call sensitivity analysis using model in class to run SAlib sobol analysis
    def sensitivity_analysis(self, sample_size):
        """Using the self variables, the sobol sensitivity analysis is used from SALib.

        These parameters are used to create the saltelli sampling problem that is then used
        in order to create a new instance of ODEModel and provides a solution.
        """
        self.problem = {
            'num_vars': self.return_string_list()[1],
            'names': self.return_string_list()[0],
            'bounds': self.param_ranges}
        # call saltelli sampling values to feed into a new model instance
        param_values = saltelli.sample(self.problem, 2**sample_size, calc_second_order = False)

        # bundle arguments needed for each model run
        args_list = [(params, self.init_cond, self.t_span, self.t_eval) for params in param_values]

        # detect available CPU cores
        num_cores = mp.cpu_count()
        print(f"...evaluating {len(param_values)} parameter sets across {num_cores} cores...")

        # spawn a pool of worker processes to evaluate models simultaneously
        with mp.Pool(processes=num_cores) as pool:
            model_out = list(
                # progress bar
                tqdm(
                    pool.imap(_solve_for_params, args_list, chunksize = 100), 
                    total=len(args_list), 
                    desc="Solving ODEs", 
                    unit="sim",
                    smoothing = 0.1
                )
            )

        # convert the results back to a numpy array for SALib
        model_out = np.array(model_out)

        sobol_analysis = sobol.analyze(self.problem, 
                                       model_out, 
                                       calc_second_order = False,
                                       n_processors = num_cores,
                                       print_to_console = True)
        return sobol_analysis


    # using values from sensitivity analysis, print out dataframe of first and total order sobol
    def sobol_dataframe_output(self, sobol_analysis):
        """Use sensitivity analysis result to give a dataframe and visualization.
        Args:
            sobol_analysis: uses sensitivity_analysis() to provide results to be used.
        Returns:
            None.
        """
        total = sobol_analysis['ST']
        conf = sobol_analysis['ST_conf']

        param_string = self.return_string_list()[0]
        df = pd.DataFrame({
            'Params': param_string,
            'Sobol First Order (S1)': sobol_analysis['S1'],
            'Sobol Total (ST)': total,
            'Sobol Total Conf': conf,
            'High Accuracy': np.where((total / conf) < 10, '--', 'Yes'),
            'Below Relative Cutoff': np.where(total < (max(total)*0.01), 'Yes', 'No')
            }).sort_values(by='Sobol Total (ST)', ascending=False)
        print("\n ...Sobol Analysis Table...")
        print(df.to_string(index=False))

        #fig, ax = plt.subplots(figsize=(10, 6))
        #df.plot(kind='bar', x='Params', y=['Sobol First Order (S1)', 'Sobol Total (ST)'], ax=ax)
        #ax.set_ylabel('Sensitivity Index Value')
        #ax.set_title('Sobol Sensitivity Analysis (Alzheimer\'s ODE Model)')
        #plt.tight_layout()
        #plt.show()

    # numerically model solution to ODE
    def solution(self):
        solution = solve_ivp(self, self.t_span, self.init_cond, t_eval=self.t_eval, method='LSODA')
        return solution

    # print results of solution
    def results(self):
        print("...calculated iterations...")
        sol = self.solution()
        print(*sol.y, sep=", ")

    def visualization(self):
        # create a grid with 2 rows and 3 columns
        fig, axs = plt.subplots(2, 3)
        # access a specific subplot using 2D indexing [row, col]
        axs[0, 0].plot(self.t_eval, self.solution().y[0])
        axs[0, 1].plot(self.t_eval, self.solution().y[1])
        axs[0, 2].plot(self.t_eval, self.solution().y[2])
        axs[1, 0].plot(self.t_eval, self.solution().y[3])
        axs[1, 1].plot(self.t_eval, self.solution().y[4])
        plt.show()
