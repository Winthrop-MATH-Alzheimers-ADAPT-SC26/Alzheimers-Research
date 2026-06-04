include("ODE.jl")
using .AlzheimerModel: build_system, make_problem, solve_model
using Statistics

function calculate_weights(; umax=(2.0, 2.0, 2.0), tspan=(0.0, 100.0))
    u1max, u2max, u3max = umax
    u1 = u1max / 2
    u2 = u2max / 2
    u3 = u3max / 2

    sys = build_system()
    prob = make_problem(sys, tspan=tspan,
        op=[sys.u₁ => u1, sys.u₂ => u2, sys.u₃ => u3])

    sol = solve_model(prob)
    _, _, _, avg_n, avg_c = mean(sol.u)

    w1 = 1.0
    I1 = u1^2 * (tspan[2] - tspan[1])
    I2 = u2^2 * (tspan[2] - tspan[1])
    I3 = u3^2 * (tspan[2] - tspan[1])
    K = w1 * avg_n

    w2 = K / avg_c
    w3 = K / I1
    w4 = K / I2
    w5 = K / I3

    return (w1=1.0, w2=w2, w3=w3, w4=w4, w5=w5)
end
