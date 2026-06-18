module Project

include("model/model.jl")
include("model/parameters.jl")

include("solvers/ODE.jl")
include("solvers/FBS-2.jl")

# include("PlottingTools.jl")
include("WeightCalculations.jl")

using .Model: sys
using .Parameters: build_param_vector
using .ODE: print_final
using .FBS: forward_backward_sweep, FBSParams
using .Weight: calculate_weights
# using .PlottingTools: ControlsPlot, SolPlotCombined, SolPlotHorizontal, make_plot

export sys, print_final,
    build_param_vector, forward_backward_sweep,
    calculate_weights, FBSParams
# ControlsPlot, SolPlotCombined, SolPlotHorizontal, make_plot

end
