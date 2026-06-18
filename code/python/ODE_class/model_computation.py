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
    temp_model = ODEModel(current_params, init_cond, t_span, t_eval)
    sol = temp_model.solution()
    return sol.y[4, -1]

# define the class for system of equations
class ODEModel:
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
        # params list is ordered: # a1, a2, b1, b2, c1, c2, c3,
        # d1, d2, e1, e2, e3, k1, k2, k3, k4, k5, sigma, r.
        self.a1 = params[0]
        self.a2 = params[1]
        self.b1 = params[2]
        self.b2 = params[3]
        self.c1 = params[4]
        self.c2 = params[5]
        self.c3 = params[6]
        self.d1 = params[7]
        self.d2 = params[8]
        self.e1 = params[9]
        self.e2 = params[10]
        self.e3 = params[11]
        self.k1 = params[12]
        self.k2 = params[13]
        self.k3 = params[14]
        self.k4 = params[15]
        self.k5 = params[16]
        self.sigma = params[17]
        self.sigma2 = params[18]
        self.r = params[19]
        #self.u1 = params[19]
        #self.u2 = params[20]
        #self.u3 = params[21]

        # for modeling solution space
        self.init_cond = initial_conditions
        self.t_span = t_span
        self.t_eval = t_eval

        self.param_ranges = ParameterEstimation(params).get_param_ranges()
        # test variables for sensitivity analysis.

    # used in sobol sensitivity analysis for problem dict
    def return_string_list(self):
        """Using parameter strings, return list of parameter names and length of list."""
        params_string = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2',
                         'c3', 'd1', 'd2', 'e1', 'e2', 'e3',
                         'k1', 'k2', 'k3', 'k4', 'k5', 'sig',
                         'sig2', 'r']#, 'u1', 'u2', 'u3'
        num_params = len(params_string)
        return [params_string, num_params]

    # use call to allow for solutions to be computed later
    def __call__(self, t, y):
        """Use this throughout to call ODE model; useful for using self keyword within class methods."""
        ABeta, Ca, Tau, N, C = y
        dABeta = self.a1 + self.a2 * (Ca / (Ca + self.sigma)) - self.k1 * ABeta #- self.u1 * ABeta
        dCa = self.b1 + self.b2 * ABeta - self.k2 * Ca  #- self.u2 * Ca
        dTau = self.c1 + self.c2 * ABeta + self.c3 * (Ca / (Ca + self.sigma2)) - self.k3 * Tau #- self.u3 * Tau
        dN = self.d1 + self.d2 * Tau - self.k4 * N
        dC = self.e1 + self.e2 * N * self.r + self.e3 * Tau - self.k5 * C
        return [dABeta, dCa, dTau, dN, dC]

    # call sensitivity analysis using model in class to run SAlib sobol analysis
    def sensitivity_analysis(self):
        """Using the self variables, the sobol sensitivity analysis is used from SALib.

        These parameters are used to create the saltelli sampling problem that is then used
        in order to create a new instance of ODEModel and provides a solution.
        """
        self.problem = {
            'num_vars': self.return_string_list()[1],
            'names': self.return_string_list()[0],
            'bounds': self.param_ranges}
        # call saltelli sampling values to feed into a new model instance
        param_values = saltelli.sample(self.problem, 2**12, calc_second_order = False)

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
                    unit="sim"
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
            'Below Relative Cutoff': np.where(total < (max(total)*0.01), '--', 'No')
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


