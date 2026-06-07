module Weight

export calculate_weights

using Statistics: mean
using DifferentialEquations: ODEProblem, solve, Tsit5
using ..Parameters: Params, build_param_vector
using ..Model: sys

function calculate_weights(x0, umax, tspan)
    u1max, u2max, u3max = umax
    u1 = u1max / 2
    u2 = u2max / 2
    u3 = u3max / 2
    p = Params(u₁=u1, u₂=u2, u₃=u3)
    pvec

    prob = ODEProblem(sys, x0, tspan, p)

    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-6)
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

    return (w1=w1, w2=w2, w3=w3, w4=w4, w5=w5)
end

end
