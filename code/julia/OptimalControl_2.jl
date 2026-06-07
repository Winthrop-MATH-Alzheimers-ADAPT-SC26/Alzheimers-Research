includet("modules/Project.jl")
using .Project

# Problem setup
const u1max = 0.12012
const u2max = 0.12012
const u3max = 0.12012
const x0 = (0.0, 100.0, 0.0, 0.0, 0.0)
const tspan = (0.0, 100.0)
const umax  = (u1max, u2max, u3max)

# Weights
weights = calculate_weights(x0, umax, tspan)

# Solve optimal control problem (FBS)
result = forward_backward_sweep(
    x0=[0.0, 100.0, 0.0, 0.0, 0.0],
    tspan=tspan,
    umax=umax,
    weights=weights,
    max_iter=500,
    tol=1e-6
)

# Plotting

fig1 = make_plot(ControlsPlot(), result)
display(fig1)

fig2 = make_plot(SolPlotCombined(), result)
display(fig2)

fig3 = make_plot(SolPlot(), result)
display(fig3)
