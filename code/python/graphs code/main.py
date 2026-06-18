import sys
import numpy as np
import matplotlib
import multiprocessing as mp 
from tqdm import tqdm

# use optimal PyPlot backend for OS 
if sys.platform == 'darwin':  # macos
    matplotlib.use('MacOSX')
elif sys.platform.startswith('linux'): # cachyos
    matplotlib.use('QtAgg')

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# systems of ODEs for solve_ivp
def ODEsystem(t, z, a1, a2, sigma1, k1, u1, b1, b2, k2, u2, c1, c2, c3, sigma2, k3, u3, d1, d2, k4, e1, e2, R, e3, k5):
    Ab, Ca, tau, N, C = z

    dAbdt = a1 + a2 * (Ca / (Ca + sigma1)) - (k1 * Ab) - (u1 * Ab)
    dCadt = b1 + (b2 * Ab) - (k2 * Ca) - (u2 * Ca)
    dtaudt = c1 + (c2 * Ab) + c3 * (Ca / (Ca + sigma2)) - (k3 * tau) - (u3 * tau)
    dNdt = d1 + (d2 * tau) - (k4 * N)
    dCdt = e1 + (e2 * N * R) + (e3 * tau) - (k5 * C)

    return [dAbdt, dCadt, dtaudt, dNdt, dCdt]

# used for treatment efficacy projections
def get_final_C(u1_val, u2_val, u3_val):
    params = (a1, a2, sigma1, k1, u1_val, b1, b2, k2, u2_val, c1, c2, c3, sigma2, k3, u3_val, d1, d2, k4, e1, e2, R, e3, k5)
    sol = solve_ivp(
        fun = ODEsystem,
        t_span = t_span,
        y0 = initials,
        args = params,
        method = 'LSODA'
    )

    # return end C from simulation
    return sol.y[4][-1]

# used for fsolving equilibrium
def find_equilibrium_sys(bestguess):
    Ab, Ca, tau, N, C = bestguess

    dAbdt = a1 + a2 * (Ca / (Ca + sigma1)) - (k1 * Ab) - (u1 * Ab)
    dCadt = b1 + (b2 * Ab) - (k2 * Ca) - (u2 * Ca)
    dtaudt = c1 + (c2 * Ab) + c3 * (Ca / (Ca + sigma2)) - (k3 * tau) - (u3 * tau)
    dNdt = d1 + (d2 * tau) - (k4 * N)
    dCdt = e1 + (e2 * N * R) + (e3 * tau) - (k5 * C)

    return [dAbdt, dCadt, dtaudt, dNdt, dCdt]

# one arg for multiprocessing
def _wrapper_get_final_C(args):
    u1, u2, u3 = args
    return get_final_C(u1, u2, u3)

# PARAMETERS ----------------

# A beta parameters
a1 = 65641 / 10000
a2 = 15778.463 / 10000
sigma1 = 100
k1 = 3035.98 / 10000
u1 = 0

# Ca parameters
b1 = 315569260
b2 = 6311385.2
k2 = 3155692.6
u2 = 0

# Tau parameters
c1 = 52.2958
c2 = 1.78367
c3 = 25
sigma2 = 100
k3 = 10.9999 / 6
u3 = 0

# N parameters
d1 = 0.07176
d2 = 0.398406
k4 = 0.3588

# C parameters
e1 = 146.308032 / 10000
e2 = 86.84 / 10000
R = 1
e3 = 199.16 / 10000
k5 = 0.8684 / 10

# ---------------------------

# initial values for Ab, Ca, tau, N, C
initials = [0, 100, 0, 0, 0]
t_span = (0, 50)
t_eval = np.linspace(t_span[0], t_span[1], 500)

if __name__ == '__main__':
    params = (a1, a2, sigma1, k1, u1, b1, b2, k2, u2, c1, c2, c3, sigma2, k3, u3, d1, d2, k4, e1, e2, R, e3, k5)

    # find equilibrium from best guess
    Ab_eq, Ca_eq, tau_eq, N_eq, C_eq = fsolve(find_equilibrium_sys, [20, 150, 60, 70, 20])
    print(f"\nEquilibrium: \nA beta: {Ab_eq:.2f} \nCalcium: {Ca_eq:.2f} \nTau: {tau_eq:.2f} \nNeuron Loss: {N_eq:.2f} \nCognitive Decline: {C_eq:.2f}\n")

    # OVERALL SIMULATION

    # run simulation forward
    print('...solving first system...')
    sol = solve_ivp(
        fun = ODEsystem,
        t_span = t_span,
        y0 = initials,
        args = params,
        t_eval = t_eval,
        method = 'LSODA'
    )

    # store results from simulation
    time_points = sol.t
    solAb = sol.y[0]
    solCa = sol.y[1]
    soltau = sol.y[2]
    solN = sol.y[3]
    solC = sol.y[4]

    plt.figure(figsize = (10, 5))
    plt.plot(time_points, solAb, label = "A beta plaques (μM)", color = '#21D2FF', lw = 2)
    plt.plot(time_points, solCa, label = "Calcium (nM)", color = '#FF2164', lw = 2)
    plt.plot(time_points, soltau, label = "Tau-p (μM)", color = '#FFBC21', lw = 2)
    plt.plot(time_points, solN, label = "Neurons Lost (Billions)", color = '#3FFF21', lw = 2)
    plt.plot(time_points, solC, label = "Cognitive Loss (MMSE score decrease)", color = '#7221FF', lw = 2)
    plt.xlabel('Years After Age 50')
    plt.ylabel('Population')
    plt.title('Alzheimer\'s ODE Simulation')
    leg = plt.legend(loc="best")
    leg.set_draggable(True)
    plt.grid(True)
    plt.show()

    # TREATMENT GRAPHS

    # ranges for the Us
    u1_range = np.linspace(0, 10, 1000)
    u2_range = np.linspace(0, 10000000, 1000)
    u3_range = np.linspace(0, 10, 1000)

    # store final C values
    C_end_u1 = []
    C_end_u2 = []
    C_end_u3 = []

    # packed args
    args1 = [(u, 0, 0) for u in u1_range]
    args2 = [(0, u, 0) for u in u2_range]
    args3 = [(0, 0, u) for u in u3_range]

    num_cores = max(1, mp.cpu_count() - 2)
    print(f'...starting treatment projections across {num_cores} cpus...')

    with mp.Pool(processes = num_cores) as pool:
        C_end_u1 = list(
            tqdm(
                pool.imap(_wrapper_get_final_C, args1),
                total = len(args1),
                desc = 'Evaluating u1',
                unit = 'sim'
            )
        )

        C_end_u2 = list(
            tqdm(
                pool.imap(_wrapper_get_final_C, args2),
                total = len(args2),
                desc = 'Evaluating u2',
                unit = 'sim'
            )
        )

        C_end_u3 = list(
            tqdm(
                pool.imap(_wrapper_get_final_C, args3),
                total = len(args3),
                desc = 'Evaluating u3',
                unit = 'sim'
            )
        )


    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # final C vs u1
    axes[0].plot(u1_range, C_end_u1, color='#21D2FF', lw=2)
    axes[0].set_ylim(0, 21)
    axes[0].set_xlabel('u1 (A beta treatment rate)')
    axes[0].set_ylabel('Final Cognitive Loss (MMSE score decrease)')
    axes[0].set_title('Impact of u1 on C at t=50')
    axes[0].grid(True)

    # final C vs u2
    axes[1].plot(u2_range, C_end_u2, color='#FF2164', lw=2)
    axes[1].set_ylim(0, 21)
    axes[1].set_xlabel('u2 (Calcium treatment rate)')
    axes[1].set_ylabel('Final Cognitive Loss (MMSE score decrease)')
    axes[1].set_title('Impact of u2 on C at t=50')
    axes[1].grid(True)

    # final C vs u3
    axes[2].plot(u3_range, C_end_u3, color='#FFBC21', lw=2)
    axes[2].set_ylim(0, 21)
    axes[2].set_xlabel('u3 (Tau treatment rate)')
    axes[2].set_ylabel('Final Cognitive Loss (MMSE score decrease)')
    axes[2].set_title('Impact of u3 on C at t=50')
    axes[2].grid(True)

    plt.tight_layout()
    plt.show()

    # PHASE PORTRAIT

    plt.figure(figsize = (8, 6))

    Ab_vals = np.linspace(-5, 35, 100)
    Ca_nullcline = (b1 + b2 * Ab_vals) / (k2 + u2)
    plt.plot(Ab_vals, Ca_nullcline, color="#5C3BFF", linestyle='--', lw=2, label='Calcium Nullcline')

    starting_points = [
        (0, 200), (10, 200), (20, 200), (30, 200),
        (0, 100), (10, 100), (20, 100), (30, 100),
        (2.5, 150), (27.5, 150)
    ]

    t_span_phase = (0, 50)
    t_eval_phase = np.linspace(t_span_phase[0], t_span_phase[1], 1000)

    print('...starting trajectory simulations...')
    for i, (Ab0, Ca0) in enumerate(starting_points):
        # initial conditions
        init_phase = [Ab0, Ca0, 0, 0, 0] 
        
        sol_phase = solve_ivp(
            fun = ODEsystem,
            t_span = t_span_phase,
            y0 = init_phase,
            args = params,
            t_eval = t_eval_phase,
            method = 'LSODA'
        )
        
        # plot the path A-beta and Calcium take
        label = 'Actual Trajectories' if i == 0 else ""
        plt.plot(sol_phase.y[0], sol_phase.y[1], color="#FF5B3B", lw=1.5, label=label)
        
        # mark where the trajectory started
        plt.scatter([Ab0], [Ca0], color = 'black', s = 20, zorder=5)

    plt.scatter([Ab_eq], [Ca_eq], color = '#3BFF5C', edgecolor = 'black', s = 25, zorder = 10, label=f'Equilibrium ({Ab_eq:.1f}, {Ca_eq:.1f})')

    plt.title('Trajectories of A beta and Calcium')
    plt.xlabel('A beta (μM)')
    plt.ylabel('Calcium (nM)')
    plt.xlim(-1, 31)
    plt.ylim(95, 205)
    leg = plt.legend(loc="best")
    leg.set_draggable(True)
    plt.grid(True, linestyle = '--', alpha = 0.6)

    plt.show()