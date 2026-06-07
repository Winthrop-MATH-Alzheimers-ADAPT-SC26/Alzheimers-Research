module PlottingTools

export ControlsPlot,
    SolPlotCombined,
    SolPlot,
    make_plot

using ..Model: sys
using CairoMakie

abstract type PlotType end

struct ControlsPlot <: PlotType end
struct SolPlotCombined <: PlotType end
struct SolPlot <: PlotType end


# Makie theme
set_theme!(Theme(
    fontsize=14,
    linewidth=2,
    Axis=(
        xgridvisible=true,
        ygridvisible=true,
        titlesize=16,
        xlabelsize=13,
    )
))

# Plot controls alongside solution
function make_plot(::ControlsPlot, result)

    t = result.t
    sol = result.sol

    Aβ = sol[sys.Aβ]
    Ca = sol[sys.Ca]
    τ = sol[sys.τ]
    N = sol[sys.N]
    C = sol[sys.C]

    u1 = result.controls.u1
    u2 = result.controls.u2
    u3 = result.controls.u3

    fig = Figure(size=(1200, 800))

    left = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()

    Axis(left[1, 1], title="Aβ")
    lines!(t, Aβ)
    Axis(left[1, 2], title="Ca")
    lines!(t, Ca)
    Axis(left[2, 1], title="τ")
    lines!(t, τ)
    Axis(left[2, 2], title="N")
    lines!(t, N)

    axC = Axis(right[1, 1], title="C")
    lines!(axC, t, C)

    axU = Axis(right[2, 1], title="Controls")
    lines!(axU, t, u1, label="u₁")
    lines!(axU, t, u2, label="u₂")
    lines!(axU, t, u3, label="u₃")

    axislegend(axU)

    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)

    return fig
end

# Single combined state plot
function make_plot(::SolPlotCombined, result)

    t = result.t
    sol = result.sol

    fig = Figure(size=(900, 600))
    ax = Axis(fig[1, 1], title="State Variables", xlabel="years")

    lines!(ax, t, sol[sys.Aβ], label="Aβ")
    lines!(ax, t, sol[sys.Ca], label="Ca")
    lines!(ax, t, sol[sys.τ], label="τ")
    lines!(ax, t, sol[sys.N], label="N")
    lines!(ax, t, sol[sys.C], label="C")

    axislegend(ax)

    return fig
end

# Individual state plots
function make_plot(::SolPlot, result)

    t = result.t
    sol = result.sol

    names = ["Aβ", "Ca", "τ", "N", "C"]
    vars = (sys.Aβ, sys.Ca, sys.τ, sys.N, sys.C)

    fig = Figure(size=(1200, 800))

    for i in 1:5
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1

        ax = Axis(fig[row, col], title=names[i], xlabel="years")
        lines!(ax, t, sol[vars[i]])
    end

    return fig
end

end
