module FBS

export forward_backward_sweep

using DifferentialEquations
using Interpolations
using LinearAlgebra
using ModelingToolkit: unknowns
using ..Model: state_rhs!, costate_rhs!

# --- Forward-Backward Sweep ─── #

function forward_backward_sweep(sys, pvec, tspan, umax, weights;
    n=365,
    max_iter=200,
    tol=1e-6
)
    t0, tf = tspan
    n = n * Int(tf - t0)
    u₁max, u₂max, u₃max = umax
    w₁, w₂, w₃, w₄, w₅ = weights
    ts = collect(range(t0, tf, length=n))

    pvec_dict = Dict(pvec)
    x0 = [pvec_dict[sys.Aβ], pvec_dict[sys.Ca], pvec_dict[sys.τ], pvec_dict[sys.N], pvec_dict[sys.C]]

    # Strip state keys from pvec for use with plain RHS functions
    unknown_syms = Set(Symbol(u) for u in unknowns(sys))
    pd = NamedTuple(Symbol(k) => v for (k, v) in pvec if Symbol(k) ∉ unknown_syms)


    # Step 1 — initialize controls to zero
    u1 = zeros(n)
    u2 = zeros(n)
    u3 = zeros(n)

    fwd_sol = nothing
    bwd_sol = nothing

    for iter in 1:max_iter

        u1_old, u2_old, u3_old = copy(u1), copy(u2), copy(u3)

        # Build interpolants from current control vectors
        u1fun = linear_interpolation(ts, u1, extrapolation_bc=Flat())
        u2fun = linear_interpolation(ts, u2, extrapolation_bc=Flat())
        u3fun = linear_interpolation(ts, u3, extrapolation_bc=Flat())


        # Step 3 — forward solve
        fwd_p = merge(pd, (u₁=u1fun, u₂=u2fun, u₃=u3fun))
        fwd_prob = ODEProblem(state_rhs!, x0, tspan, fwd_p)
        fwd_sol = solve(fwd_prob, abstol=1e-10, reltol=1e-8)


        # Step 4 — backward solve (tf → t0, λᵢ(T) = 0)
        bwd_p = merge(pd, (
            state_sol=fwd_sol,
            u₁=u1fun, u₂=u2fun, u₃=u3fun,
            w₁=w₁, w₂=w₂, w₃=w₃, w₄=w₄, w₅=w₅
        ))

        bwd_prob = ODEProblem(costate_rhs!, zeros(5), (tf, t0), bwd_p)
        bwd_sol = solve(bwd_prob, saveat=ts, abstol=1e-10, reltol=1e-8)

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

            u1_new[i] = clamp(λ1 * Aβv / (2 * w₃), 0.0, u₁max)
            u2_new[i] = clamp(λ2 * Cav / (2 * w₄), 0.0, u₂max)
            u3_new[i] = clamp(λ3 * τv / (2 * w₅), 0.0, u₃max)
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

    # Final solve with converged controls
    u1fun = linear_interpolation(ts, u1, extrapolation_bc=Flat())
    u2fun = linear_interpolation(ts, u2, extrapolation_bc=Flat())
    u3fun = linear_interpolation(ts, u3, extrapolation_bc=Flat())

    fwd_p = merge(pd, (u₁=u1fun, u₂=u2fun, u₃=u3fun))
    fwd_prob = ODEProblem(state_rhs!, x0, tspan, fwd_p)
    fwd_sol = solve(fwd_prob, saveat=ts, abstol=1e-10, reltol=1e-8)

    bwd_p = merge(pd, (
        state_sol=fwd_sol,
        u₁=u1fun, u₂=u2fun, u₃=u3fun,
        w₁=w₁, w₂=w₂, w₃=w₃, w₄=w₄, w₅=w₅
    ))
    bwd_prob = ODEProblem(costate_rhs!, zeros(5), (tf, t0), bwd_p)
    bwd_sol = solve(bwd_prob, abstol=1e-10, reltol=1e-8)

    return (
        sol=fwd_sol, t=ts, controls=(
            u1=u1,
            u2=u2,
            u3=u3
        ), costate_solution=bwd_sol
    )
end

end
