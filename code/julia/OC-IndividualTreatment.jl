include("modules/Project.jl")
includet("modules/PlottingTools.jl")
using .Project
using .PlottingTools
using CairoMakie: save

pvec = build_param_vector()
tspan = (0.0, 50.0)
umax = (10.0, 1.0e7, 10.0)
weights1 = calculate_weights(sys, pvec, tspan, (umax[1], 0.0, 0.0))
weights2 = calculate_weights(sys, pvec, tspan, (0.0, umax[2], 0.0))
weights3 = calculate_weights(sys, pvec, tspan, (0.0, 0.0, umax[3]))

p1 = FBSParams(
    u₁max=umax[1], u₂max=0.0, u₃max=0.0,
    w₁=weights1.w1, w₂=weights1.w2, w₃=weights1.w3, w₄=weights1.w4, w₅=weights1.w5,
    max_iter=1000, relax=0.5, verbose=false
)
p2 = FBSParams(
    u₁max=0.0, u₂max=umax[2], u₃max=0.0,
    w₁=weights2.w1, w₂=weights2.w2, w₃=weights2.w3, w₄=weights2.w4, w₅=weights2.w5,
    max_iter=1000, relax=0.5, verbose=false
)
p3 = FBSParams(
    u₁max=0.0, u₂max=0.0, u₃max=umax[3],
    w₁=weights3.w1, w₂=weights3.w2, w₃=weights3.w3, w₄=weights3.w4, w₅=weights3.w5,
    max_iter=1000, relax=0.5, verbose=false
)

r1 = forward_backward_sweep(sys, pvec, tspan, p1);
r2 = forward_backward_sweep(sys, pvec, tspan, p2);
r3 = forward_backward_sweep(sys, pvec, tspan, p3);

fig1 = make_plot(ControlsPlotSeperate(), r1, r2, r3)

fig2 = make_plot(ControlsPlotSeperate2(), r1, r2, r3)

fig3 = make_plot(ControlsPlotSeperate3(), r1, r2, r3)

