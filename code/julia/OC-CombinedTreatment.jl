include("modules/Project.jl")
includet("modules/PlottingTools.jl")
using .Project
using .PlottingTools

pvec = build_param_vector()
tspan = (0.0, 50.0)
umax = (10.0, 1.0e+7, 10.0)
weights = calculate_weights(sys, pvec, tspan, umax)

params = FBSParams(
    u₁max=umax[1], u₂max=umax[2], u₃max=umax[3],
    w₁=weights.w1, w₂=weights.w2, w₃=weights.w3, w₄=weights.w4, w₅=weights.w5,
    max_iter=1000, relax=0.5
)

# Solve optimal control problem (FBS)
result = forward_backward_sweep(sys, pvec, tspan, params);

# Plotting
fig1 = make_plot(ControlsPlot(), result)

fig2 = make_plot(SolPlotCombined(), result)

fig3 = make_plot(SolPlotHorizontal(), result)
