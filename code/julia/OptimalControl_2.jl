using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics
include("FBS.jl")

# State variables
@variables Aβ(t) Ca(t) τ(t) N(t) C(t)

@parameters a₁ a₂ σ k₁
@parameters b₁ b₂ k₂
@parameters c₁ c₂ c₃ k₃
@parameters d₁ d₂ k₄
@parameters e₁ e₂ e₃ k₅
@parameters R

@parameters w₁ w₂ w₃ w₄ w₅

# Control variables
@variables u₁(t) u₂(t) u₃(t)

# Costate variables
@variables λ(t)[1:5]

eqs = [
    D(Aβ) ~ a₁ + (a₂ * (Ca / (Ca + σ))) - (k₁ * Aβ) - (u₁ * Aβ),
    D(Ca) ~ b₁ + (b₂ * Aβ) - (k₂ * Ca) - (u₂ * Ca),
    D(τ) ~ c₁ + (c₂ * Aβ) + (c₃ * Ca) - (k₃ * τ) - (u₃ * τ),
    D(N) ~ d₁ + (d₂ * τ) - (k₄ * N),
    D(C) ~ e₁ + (e₂ * N * R) + (e₃ * τ) - (k₅ * C)
]

# Defining states for costate equations
states = [Aβ, Ca, τ, N, C]

# f(x,u)
rhs = [eq.rhs for eq in eqs]

# Running cost
L = w₁ * C + w₂ * N + w₃ * u₁^2 + w₄ * u₂^2 + w₅ * u₃^2

# Hamiltonian
H = L + sum(λ[i] * rhs[i] for i in 1:5)

# Costate equations
adjoint_eqs = [
    D(λ[i]) ~ -expand_derivatives(Symbolics.derivative(H, states[i]))
    for i in 1:5
]

function costate_rhs!(dλ, λ, p, t)

    Ca = p.state_sol(t; idxs=2)

    u1 = p.u1(t)
    u2 = p.u2(t)
    u3 = p.u3(t)

    λ1, λ2, λ3, λ4, λ5 = λ

    dλ[1] =
        -b₂ * λ2 -
        c₂ * λ3 +
        (k₁ + u1) * λ1

    dλ[2] =
        -c₃ * λ3 +
        (k₂ + u2) * λ2 -
        λ1 * (
            a₂ / (Ca + σ)
            -
            a₂ * Ca / (Ca + σ)^2
        )

    dλ[3] =
        -d₂ * λ4 -
        e₃ * λ5 +
        (k₃ + u3) * λ3

    dλ[4] =
        -w₂ +
        k₄ * λ4 -
        R * e₂ * λ5

    dλ[5] =
        -w₁ +
        k₅ * λ5
end

