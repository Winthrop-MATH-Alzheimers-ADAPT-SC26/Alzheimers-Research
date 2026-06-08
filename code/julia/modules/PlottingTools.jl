module PlottingTools

export ControlsPlot,
    SolPlotCombined,
    SolPlotHorizontal,
    make_plot

using CairoMakie

abstract type PlotType end

struct ControlsPlot <: PlotType end
struct SolPlotCombined <: PlotType end
struct SolPlotHorizontal <: PlotType end


# Makie theme
set_theme!(merge(Theme(
        fontsize=14,
        linewidth=2,
        Axis=(
            xgridvisible=true,
            ygridvisible=true,
            titlesize=16,
            xlabelsize=13,
        )
    ), theme_latexfonts()))

# Plot controls alongside solution
function make_plot(::ControlsPlot, result)

    t = result.t
    sol = result.sol

    Aβ = sol[1, :]
    Ca = sol[2, :]
    τ = sol[3, :]
    N = sol[4, :]
    C = sol[5, :]

    u1 = result.controls.u1
    u2 = result.controls.u2
    u3 = result.controls.u3

    fig = Figure(size=(1600, 900))

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

    fig = Figure(size=(1600, 900))
    ax = Axis(fig[1, 1], title="State Variables", xlabel="years")

    lines!(ax, t, sol[1, :], label="Aβ")
    lines!(ax, t, sol[2, :], label="Ca")
    lines!(ax, t, sol[3, :], label="τ")
    lines!(ax, t, sol[4, :], label="N")
    lines!(ax, t, sol[5, :], label="C")

    axislegend(ax)

    return fig
end

# Individual state plots
function make_plot(::SolPlotHorizontal, result)

    t = result.t
    sol = result.sol

    names = ["Aβ", "Ca", "τ", "N", "C"]

    fig = Figure(size=(1600, 900))

    for i in 1:5
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1

        ax = Axis(fig[row, col], title=names[i], xlabel="years")
        lines!(ax, t, sol[i, :])
    end

    return fig
end

end
