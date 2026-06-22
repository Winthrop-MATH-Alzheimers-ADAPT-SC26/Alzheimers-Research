using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# MTK system
@variables u(t) Aβ(t) Ca(t) τ(t) N(t) C(t)
@parameters a₁ a₂ b₁ b₂ c₁ c₂ c₃ d₁ d₂ e₁ e₂ e₃
@parameters k₁ k₂ k₃ k₄ k₅ σ₁ σ₂ R
@parameters u₁ u₂ u₃

eqs = [
    D(Aβ) ~ a₁ + a₂ * (Ca / (Ca + σ₁)) - k₁ * Aβ - u₁ * Aβ,
    D(Ca) ~ b₁ + b₂ * Aβ - k₂ * Ca - u₂ * Ca,
    D(τ) ~ c₁ + c₂ * Aβ + c₃ * (Ca / (Ca + σ₂)) - k₃ * τ - u₃ * τ,
    D(N) ~ d₁ + d₂ * τ - k₄ * N,
    D(C) ~ e₁ + e₂ * N * R + e₃ * τ - k₅ * C
]

u = [Aβ, Ca, τ, N, C]
p = [a₁, a₂, b₁, b₂, c₁, c₂, c₃, d₁, d₂, e₁, e₂, e₃,
    k₁, k₂, k₃, k₄, k₅, σ₁, σ₂, R,
    u₁, u₂, u₃]


function f(u, p) end

D(u) ~ f(u, p)

min(function L(p, n)
    sum(i -> abs(u(i, p) - ũ(i))^2 + λ * abs(u(100, p) - 1)^2, 1:n)
end)
