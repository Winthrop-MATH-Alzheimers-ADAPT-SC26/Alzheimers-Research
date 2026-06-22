using Distributed

addprocs(3; exeflags="--project")

@everywhere cd(@__DIR__)
@everywhere include("modules/SensitivityModel.jl")
@everywhere using .SensitivityModel

# IMPORTANT: also load on MAIN
include("modules/SensitivityModel.jl")
using .SensitivityModel

using GlobalSensitivity
using QuasiMonteCarlo
using ProgressMeter

# BUILD BOUNDS (main process only)
base_vals = SensitivityModel.base_vals
percent = 0.15
bounds = [(v * (1 - percent), v * (1 + percent)) for v in base_vals]

# PARALLEL EVALUATION FUNCTION
@everywhere function eval_batch(P)
    n = size(P, 2)
    Y = Matrix{Float64}(undef, 5, n)

    for j in 1:n
        Y[:, j] = collect(SensitivityModel.single_eval(@view P[:, j]))
    end

    return Y
end

# SOBOL RUN
samples = 2^14

println("Running Sobol sensitivity analysis...")

result = gsa(
    eval_batch,
    Sobol(),
    bounds;
    samples=samples,
    batch=true
)

using Printf, Statistics

function printResult(result)
    pnames = string.(SensitivityModel.params_to_test)
    onames = SensitivityModel.state_names

    println("\n=== Sobol Sensitivity Indices (sorted by total-order ST) ===\n")

    for (i, oname) in enumerate(onames)
        println("Output: $oname (max over time)")
        println(rpad("Parameter", 14), rpad("S1", 10), "ST")

        order = sortperm(result.ST[i, :], rev=true)   # most sensitive first
        for j in order
            @printf("  %-12s %-10.4f %.4f\n", pnames[j], result.S1[i, j], result.ST[i, j])
        end
        println()
    end

    println("|--- Overall Importance ---|")
    overall_importance = vec(mean(result.ST, dims=1))   # needs `using Statistics`
    for j in sortperm(overall_importance, rev=true)
        @printf("%-12s %.4f\n", pnames[j], overall_importance[j])
    end

end
