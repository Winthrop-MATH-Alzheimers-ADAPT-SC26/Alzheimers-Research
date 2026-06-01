import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. Define the ODE function: dy/dt = f(t, y)
def model(t, y):
    return -0.5 * y

# 2. Set the initial conditions and time span
y0 = [5]          # Must be structured as a list or array
t_span = (0, 10)  # (Start time, End time)
t_eval = np.linspace(0, 10, 100) # Points where you want the solution saved

# 3. Solve the ODE
sol = solve_ivp(model, t_span, y0, t_eval=t_eval)

# 4. Plot results
plt.plot(sol.t, sol.y[0], label='y(t)')
plt.xlabel('Time (t)')
plt.ylabel('y')
plt.legend()
plt.show()
