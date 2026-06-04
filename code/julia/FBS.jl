using DifferentialEquations
using Interpolations
using LinearAlgebra

# --- Parameters ─── #

const _a₁ = 65641.0
const _a₂ = 15778.463
const _b₁ = 315569260.0
const _b₂ = 6311385.2
const _c₁ = 52.2958
const _c₂ = 7.8
const _c₃ = 0.0
const _d₁ = 0.07176
const _d₂ = 0.3588
const _e₁ = 0.0
const _e₂ = 86.84
const _e₃ = 199.16
const _k₁ = 3035.98
const _k₂ = 3155692.6
const _k₃ = 10.9
const _k₄ = 0.003588
const _k₅ = 0.0
const _σ = 0.00027
const _R = 1.0

# --- State equations ─── #

function state_rhs!(dx, x, p, t)
    Aβ, Ca, τ, N, C = x

    u1 = p.u₁(t)
    u2 = p.u₂(t)
    u3 = p.u₃(t)

    dx[1] = _a₁ + (_a₂ * (Ca / (Ca + _σ))) - (_k₁ * Aβ) - (u1 * Aβ)
    dx[2] = _b₁ + (_b₂ * Aβ) - (_k₂ * Ca) - (u2 * Ca)
    dx[3] = _c₁ + (_c₂ * Aβ) + (_c₃ * Ca) - (_k₃ * τ) - (u3 * τ)
    dx[4] = _d₁ + (_d₂ * τ) - (_k₄ * N)
    dx[5] = _e₁ + (_e₂ * N * _R) + (_e₃ * τ) - (_k₅ * C)
end

# --- Costate equations ─── #
# Derived from -dλᵢ/dt = ∂H/∂xᵢ
# Terminal conditions: λᵢ(T) = 0

function costate_rhs!(dλ, λ, p, t)
    x = p.state_sol(t)
    Aβv = x[1]
    Cav = x[2]
    τv = x[3]

    u1 = p.u₁(t)
    u2 = p.u₂(t)
    u3 = p.u₃(t)
    w₁ = p.w₁
    w₂ = p.w₂
    w₃ = p.w₃
    w₄ = p.w₄
    w₅ = p.w₅

    λ1, λ2, λ3, λ4, λ5 = λ

    dλ[1] = _k₁ * λ1 + u1 * λ1 - _b₂ * λ2 - _c₂ * λ3

    dλ[2] = -(λ1 * _a₂ * _σ) / (Cav + _σ)^2 + _k₃ * λ2 + u2 * λ2 - _c₃ * λ3

    dλ[3] = _k₃ * λ3 + u3 * λ3 - _d₂ * λ4 - _e₃ * λ5

    dλ[4] = -w₁ + _k₄ * λ4 - _e₂ * _R * λ5

    dλ[5] = -w₂ + _k₅ * λ5
end

# --- Forward-Backward Sweep ─── #

function forward_backward_sweep(;
    x0=[0.01, 100.0, 0.0, 0.0, 0.0],
    tspan=(0.0, 100.0),
    umax=(2.0, 2.0, 2.0),
    weights=(1.0,1.0,1.0,1.0,1.0),
    max_iter=500,
    tol=1e-6
)
    t0, tf = tspan
    n = 365 * Int(tf)
    u₁max, u₂max, u₃max = umax
    w1, w2, w3, w4, w5 = weights
    ts = collect(range(t0, tf, length=n))

    # Step 1 — initialize controls to zero
    u1 = zeros(n)
    u2 = zeros(n)
    u3 = zeros(n)

    fwd_sol = nothing
    bwd_sol = nothing

    for iter in 1:max_iter

        u1_old = copy(u1)
        u2_old = copy(u2)
        u3_old = copy(u3)

        # Build interpolants from current control vectors
        u1fun = linear_interpolation(ts, u1, extrapolation_bc=Flat())
        u2fun = linear_interpolation(ts, u2, extrapolation_bc=Flat())
        u3fun = linear_interpolation(ts, u3, extrapolation_bc=Flat())

        # Step 3 — forward solve
        fwd_prob = ODEProblem(
            state_rhs!, x0, tspan,
            (u₁=u1fun, u₂=u2fun, u₃=u3fun)
        )
        fwd_sol = solve(fwd_prob, abstol=1e-8, reltol=1e-6)

        # Step 4 — backward solve (tf → t0, λᵢ(T) = 0)
        bwd_prob = ODEProblem(
            costate_rhs!,
            zeros(5),
            (tf, t0),
            (state_sol=fwd_sol, u₁=u1fun, u₂=u2fun, u₃=u3fun, w₁=w1, w₂=w2, w₃=w3, w₄=w4, w₅=w5)
        )
        bwd_sol = solve(bwd_prob, saveat=ts, abstol=1e-8, reltol=1e-6)

        # Step 5 — update controls using optimal characterization
        u1_new = similar(u1)
        u2_new = similar(u2)
        u3_new = similar(u3)

        for i in eachindex(ts)
            x = fwd_sol(ts[i])
            λv = bwd_sol(ts[i])   # interpolates backward sol at forward time

            Aβv = x[1]
            Cav = x[2]
            τv = x[3]

            λ1 = λv[1]
            λ2 = λv[2]
            λ3 = λv[3]

            u1_new[i] = clamp(λ1 * Aβv / (2 * _w₃), 0.0, u₁max)
            u2_new[i] = clamp(λ2 * Cav / (2 * _w₄), 0.0, u₂max)
            u3_new[i] = clamp(λ3 * τv / (2 * _w₅), 0.0, u₃max)
        end

        # Convex relaxation for stability
        u1 .= 0.5 .* u1_new .+ 0.5 .* u1_old
        u2 .= 0.5 .* u2_new .+ 0.5 .* u2_old
        u3 .= 0.5 .* u3_new .+ 0.5 .* u3_old

        # Step 6 — relative error convergence check
        err = maximum([
            norm(u1 - u1_old) / (norm(u1_old) + 1e-10),
            norm(u2 - u2_old) / (norm(u2_old) + 1e-10),
            norm(u3 - u3_old) / (norm(u3_old) + 1e-10)
        ])

        println("iter $iter  err = $err")

        if err < tol
            println("Converged at iteration $iter")
            break
        end

        if iter == max_iter
            println("Warning: did not converge in $max_iter iterations")
        end
    end

    # Final solve using converged controls

    u1fun = linear_interpolation(ts, u1, extrapolation_bc=Flat())
    u2fun = linear_interpolation(ts, u2, extrapolation_bc=Flat())
    u3fun = linear_interpolation(ts, u3, extrapolation_bc=Flat())

    fwd_prob = ODEProblem(
        state_rhs!,
        x0,
        tspan,
        (u₁=u1fun, u₂=u2fun, u₃=u3fun)
    )

    fwd_sol = solve(
        fwd_prob,
        abstol=1e-8,
        reltol=1e-6
    )

    bwd_prob = ODEProblem(
        costate_rhs!,
        zeros(5),
        (tf, t0),
        (
            state_sol=fwd_sol,
            u₁=u1fun,
            u₂=u2fun,
            u₃=u3fun, w₁=w1, w₂=w2, w₃=w3, w₄=w4, w₅=w5
        )
    )

    bwd_sol = solve(
        bwd_prob,
        abstol=1e-8,
        reltol=1e-6
    )

    return (
        t=ts,
        controls=(u1=u1, u2=u2, u3=u3),
        states=fwd_sol,
        costates=bwd_sol
    )
end
