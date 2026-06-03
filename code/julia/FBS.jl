using DifferentialEquations
using LinearAlgebra

function forward_backward_sweep(state_prob, costate_prob, optimal_control;
                                 max_iter=500, tol=1e-4)
    n    = length(state_prob.tspan[1]:1.0:state_prob.tspan[2])
    u1   = zeros(n)
    u2   = zeros(n)
    u3   = zeros(n)

    for iter in 1:max_iter
        u1_old, u2_old, u3_old = copy(u1), copy(u2), copy(u3)

        # ── Forward sweep ──────────────────────────────────────────────
        fwd_sol = solve(remake(state_prob), Rodas5P(), saveat=1.0)

        # ── Backward sweep ─────────────────────────────────────────────
        # flip tspan, solve costate backward
        bwd_sol = solve(remake(costate_prob,
                               tspan=reverse(state_prob.tspan)),
                        Rodas5P(), saveat=1.0)

        λ = reverse(bwd_sol)   # reverse to align with forward time

        # ── Update controls ────────────────────────────────────────────
        u1_new, u2_new, u3_new = optimal_control(fwd_sol, λ)

        # convex combination for stability
        u1 = 0.5u1_new + 0.5u1_old
        u2 = 0.5u2_new + 0.5u2_old
        u3 = 0.5u3_new + 0.5u3_old

        # ── Convergence ────────────────────────────────────────────────
        err = maximum([
            norm(u1 - u1_old) / (norm(u1_old) + 1e-10),
            norm(u2 - u2_old) / (norm(u2_old) + 1e-10),
            norm(u3 - u3_old) / (norm(u3_old) + 1e-10)
        ])

        println("iter $iter  err = $err")
        err < tol && (println("Converged"); break)
    end

    return u1, u2, u3
end
