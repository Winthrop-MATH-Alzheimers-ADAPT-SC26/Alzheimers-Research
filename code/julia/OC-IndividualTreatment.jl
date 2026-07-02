include("modules/Project.jl")
includet("modules/PlottingTools.jl")
using .Project
using .PlottingTools
using CairoMakie: save

pvec = build_param_vector()
tspan = (0.0, 50.0)
umax = (10.0, 1.0e7, 10.0)
weights = calculate_weights(sys, pvec, tspan, umax)

p1 = FBSParams(
    u₁max=umax[1], u₂max=0.0, u₃max=0.0,
    w₁=weights.w1, w₂=weights.w2, w₃=weights.w3, w₄=weights.w4, w₅=weights.w5,
    state_rhs=state_rhs!, costate_rhs=costate_rhs!, max_iter=10000, relax=0.5, verbose=false
)
p2 = FBSParams(
    u₁max=0.0, u₂max=umax[2], u₃max=0.0,
    w₁=weights.w1, w₂=weights.w2, w₃=weights.w3, w₄=weights.w4, w₅=weights.w5,
    state_rhs=state_rhs!, costate_rhs=costate_rhs!, max_iter=10000, relax=0.5, verbose=false
)
p3 = FBSParams(
    u₁max=0.0, u₂max=0.0, u₃max=umax[3],
    w₁=weights.w1, w₂=weights.w2, w₃=weights.w3, w₄=weights.w4, w₅=weights.w5,
    state_rhs=state_rhs!, costate_rhs=costate_rhs!, max_iter=10000, relax=0.5, verbose=false
)

r1 = forward_backward_sweep(sys, pvec, tspan, p1, n=52 * Int(tspan[end]));
r2 = forward_backward_sweep(sys, pvec, tspan, p2, n=52 * Int(tspan[end]));
r3 = forward_backward_sweep(sys, pvec, tspan, p3, n=52 * Int(tspan[end]));

# fig1 = make_plot(ControlsPlotSeperate(), r1, r2, r3)
# fig2 = make_plot(ControlsPlotSeperate2(), r1, r2, r3)
fig3 = make_plot(ConPlotSep3(), r1, r2, r3)

# --- ODE System Non-Dimensionalized
pvec_nd = build_param_vector(ND=true)
tspan_nd = Base.setindex(tspan, tspan[end] * DEFAULT_PARAMS.k₁, lastindex(tspan))
umax_nd = Tuple(x / DEFAULT_PARAMS.k₁ for x in umax)
weights_nd = calculate_weights(sys_nd, pvec_nd, tspan_nd, umax_nd)

p1_nd = FBSParams(
    u₁max=umax_nd[1], u₂max=0.0, u₃max=0.0,
    w₁=weights_nd.w1, w₂=weights_nd.w2, w₃=weights_nd.w3, w₄=weights_nd.w4, w₅=weights_nd.w5,
    state_rhs=state_rhs_nd!, costate_rhs=costate_rhs_nd!, max_iter=10000, relax=0.5, verbose=false
);
p2_nd = FBSParams(
    u₁max=0.0, u₂max=umax_nd[2], u₃max=0.0,
    w₁=weights_nd.w1, w₂=weights_nd.w2, w₃=weights_nd.w3, w₄=weights_nd.w4, w₅=weights_nd.w5,
    state_rhs=state_rhs_nd!, costate_rhs=costate_rhs_nd!, max_iter=10000, relax=0.5, verbose=false
);
p3_nd = FBSParams(
    u₁max=0.0, u₂max=0.0, u₃max=umax_nd[3],
    w₁=weights_nd.w1, w₂=weights_nd.w2, w₃=weights_nd.w3, w₄=weights_nd.w4, w₅=weights_nd.w5,
    state_rhs=state_rhs_nd!, costate_rhs=costate_rhs_nd!, max_iter=10000, relax=0.5, verbose=false
);

r1_nd = forward_backward_sweep(sys_nd, pvec_nd, tspan_nd, p1_nd, n=52 * Int(tspan[end]));
r2_nd = forward_backward_sweep(sys_nd, pvec_nd, tspan_nd, p2_nd, n=52 * Int(tspan[end]));
r3_nd = forward_backward_sweep(sys_nd, pvec_nd, tspan_nd, p3_nd, n=52 * Int(tspan[end]));

# fig1 = make_plot(ControlsPlotSeperate(), r1, r2, r3)
# fig2 = make_plot(ControlsPlotSeperate2(), r1, r2, r3)
fig3_nd = make_plot(ConPlotSep3Ver(), r1_nd, r2_nd, r3_nd; ND=true)
