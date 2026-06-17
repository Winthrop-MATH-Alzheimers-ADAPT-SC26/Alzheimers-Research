module GSA

using DifferentialEquations
using GlobalSensitivity
using DataFrames

include("../model/model.jl")
include("../model/parameters.jl")

using ..Model
using ..Parameters

export run_sobol

function run_sobol(; samples=1024)

    # Build your standard parameter vector
    p = build_param_vector()

    tspan = (0.0, 50.0)

    prob = ODEProblem(
        sys,
        [],
        tspan,
        p
    )

    #
    # Parameters we actually want to analyze
    #
    sensitivity_params = [
        sys.a₁,
        sys.a₂,
        sys.b₁,
        sys.b₂,
        sys.c₁,
        sys.c₂,
        sys.c₃,
        sys.d₁,
        sys.d₂,
        sys.e₁,
        sys.e₂,
        sys.e₃,
        sys.k₁,
        sys.k₂,
        sys.k₃,
        sys.k₄,
        sys.k₅,
        sys.σ₁,
        sys.σ₂,
        sys.R
    ]

    param_names = string.(sensitivity_params)

    #
    # Nominal values
    #
    nominal = [
        prob.ps[p]
        for p in sensitivity_params
    ]

    #
    # Bounds
    #
    lb = 0.85 .* nominal
    ub = 1.15 .* nominal

    bounds = [
        [lb[i], ub[i]]
        for i in eachindex(lb)
    ]

    #
    # Quantity of interest:
    # final value of C
    #
    function model(P)

        Y = zeros(size(P, 2))

        for i in axes(P, 2)

            local_ps = copy(prob.ps)

            for (j, param) in enumerate(sensitivity_params)
                local_ps[param] = P[j, i]
            end

            local_prob = remake(
                prob;
                p=local_ps
            )

            sol = solve(
                local_prob,
                Rodas5P(),
                abstol=1e-8,
                reltol=1e-8
            )

            if SciMLBase.successful_retcode(sol)

                Y[i] = sol[sys.C][end]

            else

                Y[i] = NaN

            end
        end

        return Y
    end

    result = gsa(
        model,
        Sobol(),
        bounds;
        samples=samples
    )

    df = DataFrame(
        Parameter=param_names,
        S1=result.S1,
        ST=result.ST
    )

    sort!(df, :ST, rev=true)

    println(df)

    return df, result
end

end
