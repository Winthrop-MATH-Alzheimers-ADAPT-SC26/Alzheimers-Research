using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics

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
