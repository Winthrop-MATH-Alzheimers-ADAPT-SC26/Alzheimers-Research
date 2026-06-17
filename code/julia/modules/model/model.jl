module Model

export sys, state_rhs!, adjoint_eqs

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D


# MTK system
@variables Aβ(t) Ca(t) τ(t) N(t) C(t)

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

@mtkcompile sys = System(eqs, t)

# --- Symbolic Solve for costate equations --- #
@variables λ(t)[1:5] w₁ w₂ w₃ w₄ w₅
states = [Aβ, Ca, τ, N, C]
rhs = [eq.rhs for eq in eqs]
L = w₁*C + w₂*N + u₁^2 + w₄*u₂^2 + w₅*u₃^2
H = L + sum(λ[i] * rhs[i] for i in 1:5)
adjoint_eqs = [
    D(λ[i]) ~ -expand_derivatives(Symbolics.derivative(H, states[i]))
    for i in 1:5
]

# State RHS
function state_rhs!(dx, x, p, t)
    Aβ, Ca, τ, N, C = x

    u1 = p.u₁(t)
    u2 = p.u₂(t)
    u3 = p.u₃(t)

    a₁, a₂ = p.a₁, p.a₂
    b₁, b₂ = p.b₁, p.b₂
    c₁, c₂, c₃ = p.c₁, p.c₂, p.c₃
    d₁, d₂ = p.d₁, p.d₂
    e₁, e₂, e₃ = p.e₁, p.e₂, p.e₃
    k₁, k₂, k₃, k₄, k₅ = p.k₁, p.k₂, p.k₃, p.k₄, p.k₅
    σ, R = p.σ, p.R

    dx[1] = a₁ + a₂ * (Ca / (Ca + σ)) - k₁ * Aβ - u1 * Aβ
    dx[2] = b₁ + b₂ * Aβ - k₂ * Ca - u2 * Ca
    dx[3] = c₁ + c₂ * Aβ + c₃ * Ca - k₃ * τ - u3 * τ
    dx[4] = d₁ + d₂ * τ - k₄ * N
    dx[5] = e₁ + e₂ * N * R + e₃ * τ - k₅ * C
end

# # Costate RHS (Hamiltonian adjoint system)
# function costate_rhs!(dλ, λ, p, t)
#     x = p.state_sol(t)
# 
#     Aβ, Ca, τ = x[1], x[2], x[3]
# 
#     u1 = p.u₁(t)
#     u2 = p.u₂(t)
#     u3 = p.u₃(t)
# 
#     a₂ = p.a₂
#     b₂, c₂, c₃ = p.b₂, p.c₂, p.c₃
#     d₂ = p.d₂
#     e₂, e₃ = p.e₂, p.e₃
#     k₁, k₂, k₃, k₄, k₅ = p.k₁, p.k₂, p.k₃, p.k₄, p.k₅
#     σ, R = p.σ, p.R
#     w₁, w₂, w₃, w₄, w₅ = p.w₁, p.w₂, p.w₃, p.w₄, p.w₅
# 
#     w₁, w₂, w₃, w₄, w₅ = p.w₁, p.w₂, p.w₃, p.w₄, p.w₅
# 
#     λ1, λ2, λ3, λ4, λ5 = λ
# 
#     dλ[1] = -b₂ * λ2 - c₂ * λ3 - (-k₁ - u₁) * λ1
#     dλ[2] = -(-k₂ - u₂) * λ2 - λ1 * (a₂ / (Ca + σ₁) + (-a₂ * Ca) / (Ca + σ₁)^2
#     dλ[3] = 
#     dλ[4] = 
#     dλ[5] = 
# end

# # Costate RHS (Hamiltonian adjoint system)
# function costate_rhs!(dλ, λ, p, t)
#     x = p.state_sol(t)
# 
#     Aβ, Ca, τ = x[1], x[2], x[3]
# 
#     u1 = p.u₁(t)
#     u2 = p.u₂(t)
#     u3 = p.u₃(t)
# 
#     a₂ = p.a₂
#     b₂, c₂, c₃ = p.b₂, p.c₂, p.c₃
#     d₂ = p.d₂
#     e₂, e₃ = p.e₂, p.e₃
#     k₁, k₂, k₃, k₄, k₅ = p.k₁, p.k₂, p.k₃, p.k₄, p.k₅
#     σ, R = p.σ, p.R
#     w₁, w₂, w₃, w₄, w₅ = p.w₁, p.w₂, p.w₃, p.w₄, p.w₅
# 
#     w₁, w₂, w₃, w₄, w₅ = p.w₁, p.w₂, p.w₃, p.w₄, p.w₅
# 
#     λ1, λ2, λ3, λ4, λ5 = λ
# 
#     dλ[1] = k₁ * λ1 + u1 * λ1 - b₂ * λ2 - c₂ * λ3
#     dλ[2] = -(a₂ * σ / (Ca + σ)^2) * λ1 + k₂ * λ2 + u2 * λ2 - c₃ * λ3
#     dλ[3] = k₃ * λ3 + u3 * λ3 - d₂ * λ4 - e₃ * λ5
#     dλ[4] = -w₁ + k₄ * λ4 - e₂ * R * λ5
#     dλ[5] = -w₂ + k₅ * λ5
# end

end
