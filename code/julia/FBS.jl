using DifferentialEquations
using Interpolations
using LinearAlgebra

function forward_backward_sweep(
    state_prob,
    costate_rhs!,
    umax;
    max_iter=500,
    tol=1e-6
)

    t0, tf = state_prob.tspan
    u₁max, u₂max, u₃max = umax
    ts = collect(range(t0, tf, length=500))
    w₃, w₄, w₅ = 1, 1, 1

    u1 = zeros(length(ts))
    u2 = zeros(length(ts))
    u3 = zeros(length(ts))

    fwd_sol = nothing
    bwd_sol = nothing

    for iter in 1:max_iter

        u1_old = copy(u1)
        u2_old = copy(u2)
        u3_old = copy(u3)

        # Current control interpolants
        u1fun = LinearInterpolation(ts, u1)
        u2fun = LinearInterpolation(ts, u2)
        u3fun = LinearInterpolation(ts, u3)

        # Forward solve
        fwd_prob = remake(
            state_prob;
            p=(
                u₁=u1fun,
                u₂=u2fun,
                u₃=u3fun
            )
        )

        fwd_sol = solve(fwd_prob, Rodas5P())

        # Transversality conditions
        λT = zeros(5)

        # Backward solve
        bwd_prob = ODEProblem(
            costate_rhs!,
            λT,
            (tf, t0),
            (
                state_sol=fwd_sol,
                u₁=u1fun,
                u₂=u2fun,
                u₃=u3fun
            )
        )

        bwd_sol = solve(bwd_prob, Rodas5P())

        # Reverse costates to match forward time
        λvals = reverse(bwd_sol.u)

        # Update controls
        u1_new = similar(u1)
        u2_new = similar(u2)
        u3_new = similar(u3)

        for i in eachindex(ts)

            x = fwd_sol(ts[i])

            Aβ = x[1]
            Ca = x[2]
            τ = x[3]

            λ1 = λvals[i][1]
            λ2 = λvals[i][2]
            λ3 = λvals[i][3]

            # Control Characteriations

            u1_new[i] = clamp(
                (λ1 * Aβ) / (2 * w₃),
                0.0,
                u₁max
            )

            u2_new[i] = clamp(
                (λ2 * Ca) / (2 * w₄),
                0.0,
                u₂max
            )

            u3_new[i] = clamp(
                (λ3 * τ) / (2 * w₅),
                0.0,
                u₃max
            )
        end

        # Relaxation
        u1 .= 0.5 .* u1_new .+ 0.5 .* u1_old
        u2 .= 0.5 .* u2_new .+ 0.5 .* u2_old
        u3 .= 0.5 .* u3_new .+ 0.5 .* u3_old

        # Convergence
        err = maximum([
            norm(u1 - u1_old) / (norm(u1_old) + 1e-10),
            norm(u2 - u2_old) / (norm(u2_old) + 1e-10),
            norm(u3 - u3_old) / (norm(u3_old) + 1e-10)
        ])

        println("iter $iter error = $err")

        if err < tol
            println("Converged")
            break
        end
    end

    return (
        controls=(u1, u2, u3),
        states=fwd_sol,
        costates=bwd_sol
    )
end
