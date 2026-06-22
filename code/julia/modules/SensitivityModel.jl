module SensitivityModel

using ModelingToolkit
using OrdinaryDiffEq
using SymbolicIndexingInterface: setp_oop

include("Project.jl")
using .Project

export single_eval, params_to_test, base_vals, state_names

# build system
const params = build_param_vector()
const tspan = (0.0, 50.0)

const prob = ODEProblem(sys, params, tspan)

# canonical MTK parameter order
const pmap = collect(ModelingToolkit.parameters(sys))

# IMPORTANT: MTK-safe parameter vector
const base_vals = [prob.ps[p] for p in pmap]

"""
We keep params_to_test ONLY for reference (not used for indexing)
"""
const params_to_test = pmap

const state_names = string.(ModelingToolkit.unknowns(sys))

# cached out-of-place setter — built once, reused on every call
const setter = setp_oop(prob, pmap)

# single evaluation
function single_eval(x::AbstractVector)
    pnew = setter(prob, collect(x))

    newprob = remake(
        prob;
        p=pnew,
        build_initializeprob=false
    )

    sol = solve(newprob)

    return (
        maximum(@view sol[1, :]),
        maximum(@view sol[2, :]),
        maximum(@view sol[3, :]),
        maximum(@view sol[4, :]),
        maximum(@view sol[5, :])
    )
end

end # module
