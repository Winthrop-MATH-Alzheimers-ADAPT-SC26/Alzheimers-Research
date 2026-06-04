import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def ODEsystem(t, z, a1, a2, sigma, k1, u1, b1, b2, k2, u2, c1, c2, c3, k3, u3, d1, d2, k4, e1, e2, R, e3, k5):
    Ab, Ca, tau, N, C = z

    dAbdt = a1 + a2 * (Ca / (Ca + sigma)) - (k1 * Ab) - (u1 * Ab)
    dCadt = b1 + (b2 * Ab) - (k2 * Ca) - (u2 * Ca)
    dtaudt = c1 + (c2 * Ab) + (c3 * Ca) - (k3 * tau) - (u3 * tau)
    dNdt = d1 + (d2 * tau) - (k4 * N)
    dCdt = e1 + (e2 * N * R) + (e3 * tau) - (k5 * C)

    return [dAbdt, dCadt, dtaudt, dNdt, dCdt]

def get_final_C(u1_val, u2_val, u3_val):
    params = (a1, a2, sigma, k1, u1_val, b1, b2, k2, u2_val, c1, c2, c3, k3, u3_val, d1, d2, k4, e1, e2, R, e3, k5)
    sol = solve_ivp(
        fun = ODEsystem,
        t_span = t_span,
        y0 = initials,
        args = params,
        method = 'LSODA'
    )

    # return end C from simulation
    return sol.y[4][-1]

def find_equilibrium(vars):
    Ab, Ca = vars
    dAb = a1 + a2 * (Ca / (Ca + sigma)) - (k1 * Ab) - (u1 * Ab)
    dCa = b1 + (b2 * Ab) - (k2 * Ca) - (u2 * Ca)
    return [dAb, dCa]

def find_equilibrium_sys(vars):
    Ab, Ca, tau, N, C = vars

    dAbdt = a1 + a2 * (Ca / (Ca + sigma)) - (k1 * Ab) - (u1 * Ab)
    dCadt = b1 + (b2 * Ab) - (k2 * Ca) - (u2 * Ca)
    dtaudt = c1 + (c2 * Ab) + (c3 * Ca) - (k3 * tau) - (u3 * tau)
    dNdt = d1 + (d2 * tau) - (k4 * N)
    dCdt = e1 + (e2 * N * R) + (e3 * tau) - (k5 * C)

    return [dAbdt, dCadt, dtaudt, dNdt, dCdt]

# A beta parameters
a1 = 65641 / 10000
a2 = 15778.463 / 10000
sigma = 100
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
c3 = 0.1
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

params = (a1, a2, sigma, k1, u1, b1, b2, k2, u2, c1, c2, c3, k3, u3, d1, d2, k4, e1, e2, R, e3, k5)

# initial values for Ab, Ca, tau, N, C
initials = [0, 100, 0, 0, 0]

Ab_eq, Ca_eq, tau_eq, N_eq, C_eq = fsolve(find_equilibrium_sys, [20, 150, 60, 70, 20])
print(f"\nEquilibrium: \nA beta: {Ab_eq:.2f} \nCalcium: {Ca_eq:.2f} \nTau: {tau_eq:.2f} \nNeuron Loss: {N_eq:.2f} \nCognitive Decline: {C_eq:.2f}\n")

t_span = (0, 50)

# OVERALL SIMULATION

t_eval = np.linspace(t_span[0], t_span[1], 500)

# run simulation forward
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
plt.xlabel('Time (Years)')
plt.ylabel('Population')
plt.title('Alzheimer\'s ODE Simulation')
plt.legend()
plt.grid(True)
plt.show()

# TREATMENT GRAPHS

# ranges for the Us
u1_range = np.linspace(0, 7.5, 500)
u2_range = np.linspace(0, 7500000, 500)
u3_range = np.linspace(0, 7.5, 500)

# store final C values
C_end_u1 = []
C_end_u2 = []
C_end_u3 = []

for u in u1_range:
    C_end_u1.append(get_final_C(u, 0, 0))

for u in u2_range:
    C_end_u2.append(get_final_C(0, u, 0))

for u in u3_range:
    C_end_u3.append(get_final_C(0, 0, u))

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# final C vs u1
axes[0].plot(u1_range, C_end_u1, color='#21D2FF', lw=2)
axes[0].set_xlabel('u1 (A beta treatment rate)')
axes[0].set_ylabel('Final Cognitive Loss (MMSE score decrease)')
axes[0].set_title('Impact of u1 on C at t=50')
axes[0].grid(True)

# final C vs u2
axes[1].plot(u2_range, C_end_u2, color='#FF2164', lw=2)
axes[1].set_xlabel('u2 (Calcium treatment rate)')
axes[1].set_ylabel('Final Cognitive Loss (MMSE score decrease)')
axes[1].set_title('Impact of u2 on C at t=50')
axes[1].grid(True)

# final C vs u3
axes[2].plot(u3_range, C_end_u3, color='#FFBC21', lw=2)
axes[2].set_xlabel('u3 (Tau treatment rate)')
axes[2].set_ylabel('Final Cognitive Loss (MMSE score decrease)')
axes[2].set_title('Impact of u3 on C at t=50')
axes[2].grid(True)

plt.tight_layout()
plt.show()

# PHASE PORTRAIT

plt.figure(figsize = (8, 6))

Ab_eq, Ca_eq = fsolve(find_equilibrium, [25, 150])

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
plt.legend()
plt.grid(True, linestyle = '--', alpha = 0.6)

plt.show()