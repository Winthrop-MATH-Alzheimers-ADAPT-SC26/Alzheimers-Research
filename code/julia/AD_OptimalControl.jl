using ModelingToolkit
using DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using InfiniteOpt, Ipopt
using Statistics
using Plots

# --- System (sys -> all symbolic) --- #

ModelingToolkit.@variables begin
    Aβ(..) = 0.0
    Ca(..) = 0.0
    τ(..) = 0.0
    N(..) = 0.0
    Cog(..) = 0.0
    u₁(..), [input = true, bounds = (0.0, 2.0)]
    u₂(..), [input = true, bounds = (0.0, 2.0)]
    u₃(..), [input = true, bounds = (0.0, 2.0)]
end

ModelingToolkit.@parameters begin
    a₁ = 65641.0
    a₂ = 15778.463
    b₁ = 315569260.0
    b₂ = 6311385.2
    c₁ = 52.2958
    c₂ = 7.8
    c₃ = 0.0
    d₁ = 0.07176
    d₂ = 0.3588
    e₁ = 0.0
    e₂ = 86.84
    e₃ = 199.16
    k₁ = 3035.98
    k₂ = 3155692.6
    k₃ = 10.9
    k₄ = 0.003588
    k₅ = 0.0
    σ = 0.00027
    R = 1.0
    w₁ = 1.0
    w₂ = 1.0
    w₃ = 1.0
end

(ts, te) = (0.0, 100000.0)

eqs = [
    D(Aβ(t)) ~ a₁ + (a₂ * (Ca(t) / (Ca(t) + σ))) - (k₁ * Aβ(t)) - (u₁(t) * Aβ(t)),
    D(Ca(t)) ~ b₁ + (b₂ * Aβ(t)) - (k₂ * τ(t)) - (u₂(t) * τ(t)),
    D(τ(t)) ~ c₁ + (c₂ * Aβ(t)) + (c₃ * Ca(t)) - (k₃ * τ(t)) - (u₃(t) * τ(t)),
    D(N(t)) ~ d₁ + (d₂ * τ(t)) - (k₄ * N(t)),
    D(Cog(t)) ~ e₁ + (e₂ * N(t) * R) + (e₃ * τ(t)) - (k₅ * Cog(t))
]

# objective: minimize cognitive decline + treatment burden
costs = [Cog(te) + w₁ * u₁(te)^2 + w₂ * u₂(te)^2 + w₃ * u₃(te)^2]

# terminal constraints — controls off at end
cons = []

@named ad_sys = System(eqs, t; costs, constraints=cons)
ad_sys = mtkcompile(ad_sys, inputs=[u₁(t), u₂(t), u₃(t)])

# --- Plain system for no-treatment comparison --- #
@named plain_sys = System(eqs, t)
plain_sys = mtkcompile(plain_sys, inputs=[u₁(t), u₂(t), u₃(t)])

# --- Initial conditions (states only, no controls) --- #
op = [
    Aβ(t) => 0.01,
    Ca(t) => 0.0,
    τ(t) => 0.0,
    N(t) => 0.0,
    Cog(t) => 0.0
]

no_treatment = solve(
    ODEProblem(plain_sys, op, (ts, te), [u₁(t) => 0.0, u₂(t) => 0.0, u₃(t) => 0.0]),
    Rodas5P(), saveat=1.0
)

# --- Build and solve optimal control problem --- #
prob = InfiniteOptDynamicOptProblem(ad_sys, op, (ts, te),
    dt=1.0,
    guesses=[
        Aβ(t) => mean(no_treatment[Aβ(t)]),
        Ca(t) => mean(no_treatment[Ca(t)]),
        τ(t) => mean(no_treatment[τ(t)]),
        N(t) => mean(no_treatment[N(t)]),
        Cog(t) => mean(no_treatment[Cog(t)]),
        u₁(t) => 0.1,
        u₂(t) => 0.1,
        u₃(t) => 0.1
    ]
)
sol = solve(prob, InfiniteOptCollocation(Ipopt.Optimizer, OrthogonalCollocation(4)),
    options=Dict(
        "nlp_scaling_method" => "gradient-based",
        "obj_scaling_factor" => 1e-6,
        "max_iter" => 5000,
        "tol" => 1e-3,
        "acceptable_tol" => 1e-2,
        "acceptable_iter" => 15,
        "print_level" => 5
    )
)
# --- No-treatment baseline --- #
# Then set them to zero in the ODE problem
no_treatment = solve(
    ODEProblem(plain_sys, op, (ts, te), [u₁(t) => 0.0, u₂(t) => 0.0, u₃(t) => 0.0]),
    Rodas5P(),
    saveat=1.0
)

# --- Plot states --- #
p1 = plot(sol.sol, idxs=[Aβ(t)], title="Aβ", xlabel="years")
p2 = plot(sol.sol, idxs=[Ca(t)], title="Ca", xlabel="years")
p3 = plot(sol.sol, idxs=[τ(t)], title="τ", xlabel="years")
p4 = plot(sol.sol, idxs=[N(t)], title="N", xlabel="years")
p5 = plot(sol.sol, idxs=[Cog(t)], title="Cognitive Decline", xlabel="years")
plot!(p5, no_treatment, idxs=[Cog(t)], label="no treatment", linestyle=:dash)

# --- Plot optimal controls --- #
p6 = plot(sol.input_sol, idxs=[u₁(t)], title="u₁ (treatment 1)", xlabel="years")
p7 = plot(sol.input_sol, idxs=[u₂(t)], title="u₂ (treatment 2)", xlabel="years")
p8 = plot(sol.input_sol, idxs=[u₃(t)], title="u₃ (treatment 3)", xlabel="years")

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(3, 3)))

println("Optimal final Cog   = ", sol.sol[Cog(t)][end])
println("Baseline final Cog  = ", no_treatment[Cog(t)][end])
