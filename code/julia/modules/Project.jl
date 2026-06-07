module Project

include("model/model.jl")
include("model/parameters.jl")

include("solvers/ODE.jl")
include("solvers/FBS.jl")

include("PlottingTools.jl")
include("WeightCalculations.jl")

using .Model: sys
using .Parameters: Params, build_param_vector
using .ODE: print_final
using .FBS: forward_backward_sweep
using .Weight: calculate_weights
using .PlottingTools: ControlsPlot, SolPlotCombined, SolPlot, make_plot

export sys, Params, print_final,
    build_param_vector, forward_backward_sweep,
    calculate_weights,
    ControlsPlot, SolPlotCombined, SolPlot, make_plot

end
