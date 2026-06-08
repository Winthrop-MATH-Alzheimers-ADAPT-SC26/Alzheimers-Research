include("modules/Project.jl")
using .Project

# Build parameter vector
const pvec = build_param_vector()

# Problem setup
u1max = 120000.0
u2max = 0.12
u3max = 0.12
tspan = (0.0, 100.0)
umax = (u1max, u2max, u3max)
#
# Weights
weights = calculate_weights(sys, pvec, tspan, umax)
#
# Solve optimal control problem (FBS)
result = forward_backward_sweep(sys, pvec, tspan, umax, weights, max_iter=1000, n=1)

# Plotting
fig1 = make_plot(ControlsPlot(), result);

fig2 = make_plot(SolPlotCombined(), result)

fig3 = make_plot(SolPlot(), result)
