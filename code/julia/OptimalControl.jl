include("modules/Project.jl")
using .Project

# Problem setup
const u1max = 0.12012
const u2max = 0.12012
const u3max = 0.12012
const tspan = (0.0, 100.0)
const umax = (u1max, u2max, u3max)

# Build parameter vector
const pvec = build_param_vector()

# Weights
weights = calculate_weights(sys, pvec, tspan, umax)

# Solve optimal control problem (FBS)
result = forward_backward_sweep(sys, pvec, tspan, umax, weights, max_iter=1000)

# Plotting
fig1 = make_plot(ControlsPlot(), result)

fig2 = make_plot(SolPlotCombined(), result)

fig3 = make_plot(SolPlot(), result)
