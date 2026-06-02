import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def ODEsystem(t, z, a1, a2, sigma, k1, u1, b1, b2, k2, u2, c1, c2, c3, k3, u3, d1, d2, k4, e1, e2, R, e3, k5):
    Ab, Ca, tau, N, C = z

    dAbdt = a1 + a2 * (Ca / (Ca + sigma)) - (k1 * Ab) - (u1 * Ab)
    dCadt = b1 + (b2 * Ab) - (k2 * Ca) - (u2 * Ca)
    dtaudt = c1 + (c2 * Ab) + (c3 * Ca) - (k3 * tau) - (u3 * tau)
    dNdt = d1 + (d2 * tau) - (k4 * N)
    dCdt = e1 + (e2 * N * R) + (e3 * tau) - (k5 * C)

    return [dAbdt, dCadt, dtaudt, dNdt, dCdt]

a1 = 65641 / 10000
a2 = 15778.463 / 10000
sigma = 0.00027
k1 = 3035.98 / 10000
u1 = 0
b1 = 315569260
b2 = 6311385.2
k2 = 3155692.6
u2 = 0
c1 = 52.2958
c2 = 1.78367
c3 = 0.1
k3 = 10.9999 / 6
u3 = 0
d1 = 0.07176
d2 = 0.398406
k4 = 0.3588
e1 = 146.308032 / 10000
e2 = 86.84 / 10000
R = 1
e3 = 199.16 / 10000
k5 = 0.8684 / 10

params = (a1, a2, sigma, k1, u1, b1, b2, k2, u2, c1, c2, c3, k3, u3, d1, d2, k4, e1, e2, R, e3, k5)

initials = [0, 100, 0, 0, 0]

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
plt.plot(time_points, solAb, label = "A beta plaques (μM)", color = 'blue', lw = 2)
plt.plot(time_points, solCa, label = "Calcium (nM)", color = 'red', lw = 2)
plt.plot(time_points, soltau, label = "Tau-p (μM)", color = 'green', lw = 2)
plt.plot(time_points, solN, label = "Neurons Lost (Billions)", color = 'orange', lw = 2)
plt.plot(time_points, solC, label = "Cognitive Loss (MMSE score)", color = 'purple', lw = 2)
plt.xlabel('Time (years)')
plt.ylabel('Population')
plt.title('Alzheimer\'s ODE Simulation')
plt.legend()
plt.grid(True)
plt.show()

# TREATMENT GRAPHS

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
axes[0].plot(u1_range, C_end_u1, color='blue', lw=2)
axes[0].set_xlabel('u1 (A beta treatment rate)')
axes[0].set_ylabel('Final Cognitive Loss (MMSE score)')
axes[0].set_title('Impact of u1 on C at t=50')
axes[0].grid(True)

# final C vs u2
axes[1].plot(u2_range, C_end_u2, color='red', lw=2)
axes[1].set_xlabel('u2 (Calcium treatment rate)')
axes[1].set_ylabel('Final Cognitive Loss (MMSE score)')
axes[1].set_title('Impact of u2 on C at t=50')
axes[1].grid(True)

# final C vs u3
axes[2].plot(u3_range, C_end_u3, color='green', lw=2)
axes[2].set_xlabel('u3 (Tau treatment rate)')
axes[2].set_ylabel('Final Cognitive Loss (MMSE score)')
axes[2].set_title('Impact of u3 on C at t=50')
axes[2].grid(True)

plt.tight_layout()
plt.show()