module FBS

export forward_backward_sweep, FBSParams

using DifferentialEquations
using Interpolations
using LinearAlgebra
using ModelingToolkit: unknowns
using ..Model: state_rhs!, costate_rhs!

# --- Solver / control tuning parameters ─────────────────────────────── #
#
# Bundles the knobs that travel together conceptually (control bounds,
# objective weights, convergence settings) and are independent of the
# model (`sys`), parameters (`pvec`), and time window (`tspan`).
#
# All fields are required except tol/max_iter/relax/verbose, which have
# sensible defaults. Construct with keywords, e.g.:
#
#   params = FBSParams(u₁max=1.0, u₂max=1.0, u₃max=1.0,
#                       w₁=1.0, w₂=1.0, w₃=1.0, w₄=1.0, w₅=1.0)
#
Base.@kwdef struct FBSParams
    u₁max::Float64
    u₂max::Float64
    u₃max::Float64
    w₁::Float64
    w₂::Float64
    w₃::Float64
    w₄::Float64
    w₅::Float64
    tol::Float64 = 1e-6
    max_iter::Int = 200
    relax::Float64 = 0.5      # weight on the NEW control in convex combination;
    # relax=1.0 recovers unrelaxed FBS (often unstable)
    verbose::Bool = true
end

# --- Internal helper: one forward + backward solve at fixed controls ── #
#
# Factored out so the in-loop solve and the final post-convergence solve
# are guaranteed to use identical logic — no duplicated/divergent copies.
#
# `bwd_saveat` lets the caller control the backward grid explicitly so the
# returned costate solution can be made consistent with `ts` (or left as a
# dense continuous interpolant when `nothing`).
function _solve_states_and_costates(
    state_rhs!, costate_rhs!,
    x0, tspan, pd,
    u1fun, u2fun, u3fun,
    weights, ts;
    bwd_saveat=nothing,
    abstol=1e-10, reltol=1e-8
)
    t0, tf = tspan
    w₁, w₂, w₃, w₄, w₅ = weights

    fwd_p = merge(pd, (u₁=u1fun, u₂=u2fun, u₃=u3fun))
    fwd_prob = ODEProblem(state_rhs!, x0, tspan, fwd_p)
    fwd_sol = solve(fwd_prob, saveat=ts, abstol=abstol, reltol=reltol)

    if fwd_sol.retcode != ReturnCode.Success
        error("Forward state solve failed with retcode=$(fwd_sol.retcode). " *
              "Costate/control update would be unreliable downstream.")
    end

    bwd_p = merge(pd, (
        state_sol=fwd_sol,
        u₁=u1fun, u₂=u2fun, u₃=u3fun,
        w₁=w₁, w₂=w₂, w₃=w₃, w₄=w₄, w₅=w₅
    ))

    bwd_prob = ODEProblem(costate_rhs!, zeros(5), (tf, t0), bwd_p)
    bwd_sol = bwd_saveat === nothing ?
              solve(bwd_prob, abstol=abstol, reltol=reltol) :
              solve(bwd_prob, saveat=bwd_saveat, abstol=abstol, reltol=reltol)

    if bwd_sol.retcode != ReturnCode.Success
        error("Backward costate solve failed with retcode=$(bwd_sol.retcode). " *
              "Check sign convention in costate_rhs! (t runs t0→tf even though " *
              "integration direction is tf→t0).")
    end

    return fwd_sol, bwd_sol
end

# --- Internal helper: optimality-condition control update ──────────── #
#
# All three controls follow the same λᵢ·xᵢ / (2wᵢ) stationarity form;
# looping over the triples makes that structural symmetry visible instead
# of repeating near-identical lines three times.
function _update_controls(fwd_sol, bwd_sol, ts, weights_uw, umax)
    u_new = ntuple(_ -> similar(ts), 3)

    for i in eachindex(ts)
        x = fwd_sol(ts[i])
        λv = bwd_sol(ts[i])
        for k in 1:3
            λk = λv[k]
            xk = x[k]
            wk = weights_uw[k]
            umaxk = umax[k]
            u_new[k][i] = (wk == 0.0) ? 0.0 : clamp(λk * xk / (2 * wk), 0.0, umaxk)
        end
    end

    # Boundary convention: control held at zero at t0. This is a modeling
    # assumption (treatment cannot start instantaneously at the initial
    # time), not a numerical patch — documented here so it isn't mistaken
    # for leftover indexing cleanup.
    for k in 1:3
        u_new[k][1] = 0.0
    end

    return u_new
end

# --- Forward-Backward Sweep ──────────────────────────────────────────── #

function forward_backward_sweep(sys, pvec, tspan, params::FBSParams; n=365)
    t0, tf = tspan
    n = n * Int(tf - t0)   # n is points-per-unit-time (e.g. daily resolution);
    # total grid size scales with the horizon length
    ts = collect(range(t0, tf, length=n))

    umax = (params.u₁max, params.u₂max, params.u₃max)
    weights = (params.w₁, params.w₂, params.w₃, params.w₄, params.w₅)
    weights_uw = (params.w₃, params.w₄, params.w₅)  # weights entering the control update

    pvec_dict = Dict(pvec)
    x0 = [pvec_dict[sys.Aβ], pvec_dict[sys.Ca], pvec_dict[sys.τ], pvec_dict[sys.N], pvec_dict[sys.C]]

    # Strip state keys from pvec for use with plain RHS functions
    unknown_syms = Set(Symbol(u) for u in unknowns(sys))
    pd = NamedTuple(Symbol(k) => v for (k, v) in pvec if Symbol(k) ∉ unknown_syms)

    # Baseline (no treatment)
    zero_p = merge(pd, (u₁=t -> 0.0, u₂=t -> 0.0, u₃=t -> 0.0))
    baseline_prob = ODEProblem(state_rhs!, x0, tspan, zero_p)
    baseline_sol = solve(baseline_prob, saveat=ts, abstol=1e-10, reltol=1e-8)

    # Step 1 — initialize controls to zero
    u1 = zeros(n)
    u2 = zeros(n)
    u3 = zeros(n)

    fwd_sol = nothing
    bwd_sol = nothing
    err_history = Float64[]   # exposed so callers can diagnose non-convergence
    # (oscillating vs. plateauing vs. diverging)

    for iter in 1:params.max_iter

        u1_old, u2_old, u3_old = copy(u1), copy(u2), copy(u3)

        u1fun = linear_interpolation(ts, u1, extrapolation_bc=Flat())
        u2fun = linear_interpolation(ts, u2, extrapolation_bc=Flat())
        u3fun = linear_interpolation(ts, u3, extrapolation_bc=Flat())

        fwd_sol, bwd_sol = _solve_states_and_costates(
            state_rhs!, costate_rhs!,
            x0, tspan, pd,
            u1fun, u2fun, u3fun,
            weights, ts;
            bwd_saveat=reverse(ts)
        )

        # Step 5 — update controls using optimal characterization
        u1_new, u2_new, u3_new = _update_controls(fwd_sol, bwd_sol, ts, weights_uw, umax)

        # Step 6 — relative error convergence check.
        err = maximum([
            norm(u1_new - u1_old) / (norm(u1_old) + 1e-10),
            norm(u2_new - u2_old) / (norm(u2_old) + 1e-10),
            norm(u3_new - u3_old) / (norm(u3_old) + 1e-10)
        ])
        push!(err_history, err)

        # Convex relaxation for stability — applied AFTER the convergence
        r = params.relax
        u1 .= r .* u1_new .+ (1 - r) .* u1_old
        u2 .= r .* u2_new .+ (1 - r) .* u2_old
        u3 .= r .* u3_new .+ (1 - r) .* u3_old

        params.verbose && println("iter $iter  err = $err")

        if err < params.tol
            println("Converged at iteration $iter")
            break
        end

        if iter == params.max_iter
            @warn "Did not converge in $(params.max_iter) iterations (final err=$err). " *
                  "Consider lowering `relax` or inspecting `err_history` for oscillation."
        end
    end

    # Final solve with converged controls — reuses the exact same helper
    # as the loop, so this can never silently drift from in-loop behavior.
    u1fun = linear_interpolation(ts, u1, extrapolation_bc=Flat())
    u2fun = linear_interpolation(ts, u2, extrapolation_bc=Flat())
    u3fun = linear_interpolation(ts, u3, extrapolation_bc=Flat())

    fwd_sol, bwd_sol = _solve_states_and_costates(
        state_rhs!, costate_rhs!,
        x0, tspan, pd,
        u1fun, u2fun, u3fun,
        weights, ts;
        bwd_saveat=reverse(ts)   # kept consistent with `ts`, unlike the original
    )

    return (
        sol=fwd_sol,
        baseline_sol=baseline_sol,
        t=ts,
        controls=(u1=u1, u2=u2, u3=u3),
        costate_solution=bwd_sol,
        err_history=err_history
    )
end

end
