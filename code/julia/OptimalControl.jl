using ModelingToolkit
using DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using InfiniteOpt, Ipopt
using Statistics
using Plots

# --- System (sys -> all symbolic) --- #

ModelingToolkit.@variables begin
    Aβ(..) = 0.0
    Ca(..) = 100.0
    τ(..) = 0.0
    N(..) = 0.0
    C(..) = 0.0
    u₁(..), [input = true, bounds = (0.0, 1.0)]
    u₂(..), [input = true, bounds = (0.0, 1.0)]
    u₃(..), [input = true, bounds = (0.0, 1.0)]
end

ModelingToolkit.@parameters begin
    a₁ = 4.380
    a₂ = 1.578
    b₁ = 315569260.0
    b₂ = 631139.0
    c₁ = 52.2958
    c₂ = 7.8
    c₃ = 0.1
    d₁ = 0.07176
    d₂ = 0.3588
    e₁ = 0.0146308032
    e₂ = 0.008684
    e₃ = 0.019916
    k₁ = 4.38
    k₂ = 3155692.6
    k₃ = 1.833316667
    k₄ = 0.3588
    k₅ = 0.08684
    σ = 0.00027
    R = 1.0
end

(ts, te) = (0.0, 50.0)

eqs = [
    D(Aβ(t)) ~ a₁ + (a₂ * (Ca(t) / (Ca(t) + σ))) - (k₁ * Aβ(t)) - (u₁(t) * Aβ(t)),
    D(Ca(t)) ~ b₁ + (b₂ * Aβ(t)) - (k₂ * Ca(t)) - (u₂(t) * Ca(t)),
    D(τ(t)) ~ c₁ + (c₂ * Aβ(t)) + (c₃ * Ca(t)) - (k₃ * τ(t)) - (u₃(t) * τ(t)),
    D(N(t)) ~ d₁ + (d₂ * τ(t)) - (k₄ * N(t)),
    D(C(t)) ~ e₁ + (e₂ * N(t) * R) + (e₃ * τ(t)) - (k₅ * C(t))
]

# objective: minimize cognitive decline + treatment burden
costs = [C(te), N(te)]

# terminal constraints — controls off at end
cons = []

@named ad_sys = System(eqs, t; costs, constraints=cons)
ad_sys = mtkcompile(ad_sys, inputs=[u₁(t), u₂(t), u₃(t)])

# --- Plain system for no-treatment comparison --- #
@named plain_sys = System(eqs, t)
plain_sys = mtkcompile(plain_sys, inputs=[u₁(t), u₂(t), u₃(t)])

# --- Initial conditions (states only, no controls) --- #
op = [
    Aβ(t) => 0.00,
    Ca(t) => 100.0,
    τ(t) => 0.0,
    N(t) => 0.0,
    C(t) => 0.0
]

no_treatment = solve(ODEProblem(plain_sys, op, (ts, te), [u₁(t) => 0.0, u₂(t) => 0.0, u₃(t) => 0.0]))

# --- Build and solve optimal control problem --- #
prob = InfiniteOptDynamicOptProblem(ad_sys, op, (ts, te),
    dt=0.1,
    guesses=[
        Aβ(t) => mean(no_treatment[Aβ(t)]),
        Ca(t) => mean(no_treatment[Ca(t)]),
        τ(t) => mean(no_treatment[τ(t)]),
        N(t) => mean(no_treatment[N(t)]),
        C(t) => mean(no_treatment[C(t)]),
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

# --- Plot states --- #
p1 = plot(sol.sol, idxs=[Aβ(t)], title="Aβ", xlabel="years")
p2 = plot(sol.sol, idxs=[Ca(t)], title="Ca", xlabel="years")
p3 = plot(sol.sol, idxs=[τ(t)], title="τ", xlabel="years")
p4 = plot(sol.sol, idxs=[N(t)], title="N", xlabel="years")
p5 = plot(sol.sol, idxs=[C(t)], title="C", xlabel="years")
plot!(p5, no_treatment, idxs=[C(t)], label="no treatment", linestyle=:dash)
# --- Plot optimal controls --- #

t_ctrl = sol.input_sol.t
u1_vals = [u[1] for u in sol.input_sol.u]
u2_vals = [u[2] for u in sol.input_sol.u]
u3_vals = [u[3] for u in sol.input_sol.u]

p6 = plot(t_ctrl, u1_vals, title="u₁ (treatment 1)", xlabel="years", ylim=(0, 1), label="u₁")
p7 = plot(t_ctrl, u2_vals, title="u₂ (treatment 2)", xlabel="years", ylim=(0, 1), label="u₂")
p8 = plot(t_ctrl, u3_vals, title="u₃ (treatment 3)", xlabel="years", ylim=(0, 1), label="u₃")

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(2, 4), size=(2496, 1664), dpi=150, plot_titlefontsize=14))
