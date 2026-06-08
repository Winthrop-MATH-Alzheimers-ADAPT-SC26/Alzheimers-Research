include("modules/Project.jl")
using .Project

# Build parameter vector
a₁ = 0.000002
a₂ = 1.5
b₁ = 5.0e-3
b₂ = 10.2
c₁ = 0.01
c₂ = 1.0
c₃ = 0.2
d₁ = 0.001
d₂ = 0.4
e₁ = 0.0
e₂ = 0.2
e₃ = 0.15
k₁ = 0.004
k₂ = 0.05
k₃ = 0.002
k₄ = 0.015
k₅ = 0.008
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
u1max = 0.0
u2max = 0.0
u3max = 0.0
tspan = (0.0, 100.0)
umax = (u1max, u2max, u3max)
#
# Weights
weights = calculate_weights(sys, pvec, tspan, umax)
#
# Test Weights
weights_new = (; weights..., w3=1.0e5, w4=1.0e-6, w5=1.0e-3)
#
# Solve optimal control problem (FBS)
result = forward_backward_sweep(sys, pvec, tspan, umax, weights_new, max_iter=1000, n=1)
#
# Plotting
fig3 = make_plot(SolPlotHorizontal(), result)

fig1 = make_plot(ControlsPlot(), result)

fig2 = make_plot(SolPlotCombined(), result)

