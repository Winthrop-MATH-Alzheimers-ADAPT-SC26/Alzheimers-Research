module Project

include("model/model.jl")
include("model/parameters.jl")

include("solvers/ODE.jl")
include("solvers/FBS.jl")

# include("PlottingTools.jl")
include("WeightCalculations.jl")

using .Model: sys, sys_nd, state_rhs!, costate_rhs!, state_rhs_nd!, costate_rhs_nd!
using .Parameters: build_param_vector, DEFAULT_PARAMS
using .ODE: print_final
using .FBS: forward_backward_sweep, FBSParams
using .Weight: calculate_weights

export sys, sys_nd, state_rhs!, costate_rhs!, state_rhs_nd!, costate_rhs_nd!,
    print_final, build_param_vector, forward_backward_sweep,
    calculate_weights, FBSParams, DEFAULT_PARAMS

println("Project Module Loaded")

end
