include("modules/Project.jl")
includet("modules/PlottingTools.jl")
using .Project
using .PlottingTools

# Build parameter vector
a₁ = 0.000002
a₂ = 1.5
b₁ = 5.0e-3
b₂ = 1.2
c₁ = 0.001
c₂ = 1.0
c₃ = 0.02
d₁ = 0.001
d₂ = 0.04
e₁ = 0.0
e₂ = 0.002
e₃ = 0.00015
k₁ = 0.004
k₂ = 0.05
k₃ = 0.002
k₄ = 0.015
k₅ = e₁ / 30
u₁ = 0.0
u₂ = 0.0
u₃ = 0.0
σ = 5.0
Aβ = 0.01
Ca = 0.05
τ = 0.01
p_kwargs = (; a₁, a₂, k₁, u₁, b₁, b₂, k₂, u₂,
    c₁, c₂, c₃, k₃, u₃, d₁, d₂, k₄,
    e₁, e₂, e₃, k₅)
#
pvec = build_param_vector(; p_kwargs...)
# Problem setup
# u1max = 2.0
# u2max = 3.0
# u3max = 2.5
tspan = (0.0, 50.0)
# umax = (u1max, u2max, u3max)
#
function get_result(sys, pvec, tspan, umax; weights_new=nothing)
    if !isnothing(weights_new)
        w3, w4, w5 = weights_new
        weights = calculate_weights(sys, pvec, tspan, umax)
        weights_new_combined = (; weights..., w3=w3, w4=w4, w5=w5)
        result = forward_backward_sweep(sys, pvec, tspan, umax, weights_new_combined, max_iter=1000, n=4)
    else
        weights = calculate_weights(sys, pvec, tspan, umax)
        result = forward_backward_sweep(sys, pvec, tspan, umax, weights, max_iter=1000, n=4)
    end
    return result
end

umax1 = (10.0, 0.0, 0.0)
umax2 = (0.0, 10.0, 0.0)
umax3 = (0.0, 0.0, 10.0)

pvec_default = build_param_vector()
r1 = get_result(sys, pvec_default, tspan, umax1)
r2 = get_result(sys, pvec_default, tspan, umax2)
r3 = get_result(sys, pvec_default, tspan, umax3)

fig4 = make_plot(ControlsPlotSeperate2(), r1, r2, r3)

umax = (10.0, 1.0e+7, 10.0)
result = get_result(sys, pvec_default, tspan, umax)

# Plotting
fig1 = make_plot(ControlsPlotSeperate(), r1, r2, r3)

fig2 = make_plot(SolPlotCombined(), result)

fig3 = make_plot(SolPlotHorizontal(), result)
