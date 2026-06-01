module AlzheimerModel

export build_system, make_problem, solve_model, plot_solution, print_final

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using Plots

# --- System (sys -> all symbolic) --- #

function build_system()
    @variables begin
        Aβ(t) = 0.01
        Ca(t) = 0.0
        τ(t) = 0.0
        N(t) = 0.0
        C(t) = 0.0
    end

    @parameters begin
        a₁ = 65641
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
        u₁ = 0.0
        u₂ = 0.0
        u₃ = 0.0
        σ = 0.00027
        R = 1.0
    end

    eqs = [
        D(Aβ) ~ a₁ + (a₂ * (Ca / (Ca + σ))) - (k₁ * Aβ) - (u₁ * Aβ),
        D(Ca) ~ b₁ + (b₂ * Aβ) - (k₂ * τ) - (u₂ * τ),
        D(τ) ~ c₁ + (c₂ * Aβ) + (c₃ * Ca) - (k₃ * τ) - (u₃ * τ),
        D(N) ~ d₁ + (d₂ * τ) - (k₄ * N),
        D(C) ~ e₁ + (e₂ * N * R) + (e₃ * τ) - (k₅ * C)
    ]

    @mtkcompile sys = ODESystem(eqs, t)
    return sys
end

# --- Problem (prob -> concrete numbers plugged in, ready to solve) --- #

function make_problem(sys; tspan=(0.0, 100.0), u₀=[], p=[])
    return ODEProblem(sys, u₀, tspan, p)
end

# --- Solve (sol -> final result, queryable across time) --- #

function solve_model(prob; solver=nothing, abstol=1e-8, reltol=1e-6)
    if solver === nothing
        return solve(prob, abstol=abstol, reltol=reltol)
    else
        return solve(prob, solver, abstol=abstol, reltol=reltol)
    end
end

# --- Output --- #

function plot_solution(sol, sys)
    p1 = plot(sol, idxs=[sys.Aβ], title="Aβ")
    p2 = plot(sol, idxs=[sys.Ca], title="Ca")
    p3 = plot(sol, idxs=[sys.τ], title="τ")
    p4 = plot(sol, idxs=[sys.N], title="N")
    p5 = plot(sol, idxs=[sys.C], title="C")
    return plot(p1, p2, p3, p4, p5, layout=(2, 3))
end

function print_final(sol, sys)
    println("Fianl values at t = ", sol.t[end])
    println("  Aβ = ", sol[sys.Aβ][end])
    println("  Ca = ", sol[sys.Ca][end])
    println("  τ  = ", sol[sys.τ][end])
    println("  N  = ", sol[sys.N][end])
    println("  C  = ", sol[sys.C][end])
end

end
