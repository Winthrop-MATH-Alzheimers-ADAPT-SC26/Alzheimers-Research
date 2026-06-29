module Weight

export calculate_weights

using Statistics: mean
using DifferentialEquations
using ModelingToolkit

function calculate_weights(sys::ModelingToolkitBase.System, pvec::Vector, tspan::NTuple{2,Float64}, umax::NTuple{3,Float64})
    u1max, u2max, u3max = umax
    u1 = u1max / 2
    u2 = u2max / 2
    u3 = u3max / 2
    t0, tf = tspan
    n = 365 * Int(tf - t0)
    ts = collect(range(t0, tf, length=n))

    pvec_u = [k => v for (k, v) in merge(Dict(pvec), Dict(sys.u₁ => u1, sys.u₂ => u2, sys.u₃ => u3))]

    prob = ODEProblem(sys, pvec_u, tspan)

    sol = solve(prob, abstol=1e-8, reltol=1e-6, saveat=ts)
    avg_n = mean(sol[sys.N])
    avg_c = mean(sol[sys.C])

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
