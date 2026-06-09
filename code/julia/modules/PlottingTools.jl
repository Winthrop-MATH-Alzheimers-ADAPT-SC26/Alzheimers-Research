module PlottingTools

export ControlsPlot,
    ControlsPlotSeperate,
    ControlsPlotSeperate2,
    SolPlotCombined,
    SolPlotHorizontal,
    make_plot

using CairoMakie

abstract type PlotType end

struct ControlsPlot <: PlotType end
struct ControlsPlotSeperate <: PlotType end
struct ControlsPlotSeperate2 <: PlotType end
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

function make_plot(::ControlsPlotSeperate, r1, r2, r3)
    t = r1.t .+ 50.0
    sol1, sol2, sol3 = r1.sol, r2.sol, r3.sol

    C1 = sol1[5, :]
    C2 = sol2[5, :]
    C3 = sol3[5, :]

    C1_base = r1.baseline_sol[5, :]
    C2_base = r2.baseline_sol[5, :]
    C3_base = r3.baseline_sol[5, :]

    u1 = r1.controls.u1
    u2 = r2.controls.u2
    u3 = r3.controls.u3

    fig = Figure(size=(1600, 900))
    left = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()

    axc1 = Axis(left[1, 1], title="C with only u₁ Treatment", xticks=50:10:100)
    lines!(axc1, t, C1, label="With Treatment")
    lines!(axc1, t, C1_base, linestyle=:dash, label="No Treatment")
    axislegend(axc1)

    axc2 = Axis(left[2, 1], title="C with only u₂ Treatment", xticks=50:10:100)
    lines!(axc2, t, C2, label="With Treatment")
    lines!(axc2, t, C2_base, linestyle=:dash, label="No Treatment")
    axislegend(axc2)

    axc3 = Axis(left[3, 1], title="C with only u₃ Treatment", xticks=50:10:100)
    lines!(axc3, t, C3, label="With Treatment")
    lines!(axc3, t, C3_base, linestyle=:dash, label="No Treatment")
    axislegend(axc3)

    axc4l = Axis(right[1, 1], title="Each Treatment", xticks=50:10:100)
    axc4r = Axis(right[1, 1], xticks=50:10:100, yaxisposition=:right)
    l1 = lines!(axc4l, t, u1, label="u₁")
    l3 = lines!(axc4l, t, u3, label="u₃")
    l2 = lines!(axc4r, t, u2, label="u₂")

    Legend(
        right[1, 1],
        [l1, l2, l3],
        ["u₁", "u₂", "u₃"],
        tellwidth=false,
        tellheight=false,
        halign=:right,
        valign=:top
    )

    return fig
end

function make_plot(::ControlsPlotSeperate2, r1, r2, r3)
    t = r1.t .+ 50.0
    sol1, sol2, sol3 = r1.sol, r2.sol, r3.sol

    C1 = sol1[5, :]
    C2 = sol2[5, :]
    C3 = sol3[5, :]

    C1_base = r1.baseline_sol[5, :]
    C2_base = r2.baseline_sol[5, :]
    C3_base = r3.baseline_sol[5, :]

    u1 = r1.controls.u1
    u2 = r2.controls.u2
    u3 = r3.controls.u3

    fig = Figure(size=(1600, 900))
    left = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()

    axc1 = Axis(left[1, 1], title="C with only u₁ Treatment", xticks=50:10:100)
    lines!(axc1, t, C1, label="With Treatment")
    lines!(axc1, t, C1_base, linestyle=:dash, label="No Treatment")
    axislegend(axc1)

    axc2 = Axis(left[2, 1], title="C with only u₂ Treatment", xticks=50:10:100)
    lines!(axc2, t, C2, label="With Treatment")
    lines!(axc2, t, C2_base, linestyle=:dash, label="No Treatment")
    axislegend(axc2)

    axc3 = Axis(left[3, 1], title="C with only u₃ Treatment", xticks=50:10:100)
    lines!(axc3, t, C3, label="With Treatment")
    lines!(axc3, t, C3_base, linestyle=:dash, label="No Treatment")
    axislegend(axc3)

    axc4 = Axis(right[1, 1], title="u₁ Treatment", xticks=50:10:100)
    lines!(axc4, t, u1, label="u₁")
    axislegend(axc4)

    axc5 = Axis(right[2, 1], title="u₂ Treatment", xticks=50:10:100)
    lines!(axc5, t, u2, label="u₂")
    axislegend(axc5)

    axc6 = Axis(right[3, 1], title="u₃ Treatment", xticks=50:10:100)
    lines!(axc6, t, u3, label="u₃")
    axislegend(axc6)

    return fig
end
# Plot controls alongside solution
function make_plot(::ControlsPlot, result)

    t = result.t .+ 50.0
    sol = result.sol

    Aβ = sol[1, :]
    Ca = sol[2, :]
    τ = sol[3, :]
    N = sol[4, :]
    C = sol[5, :]
    C_no_treatment = result.baseline_sol[5, :]

    u1 = result.controls.u1
    u2 = result.controls.u2
    u3 = result.controls.u3

    fig = Figure(size=(1600, 900))

    left = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()

    Axis(left[1, 1], title="Aβ Levels", xticks=50:10:100)
    lines!(t, Aβ)
    Axis(left[1, 2], title="Ca Levels", xticks=50:10:100)
    lines!(t, Ca)
    Axis(left[2, 1], title="τ Levels", xticks=50:10:100)
    lines!(t, τ)
    Axis(left[2, 2], title="N Levels", xticks=50:10:100)
    lines!(t, N)

    axC = Axis(right[1, 1], title="C Levels", xticks=50:10:100)
    lines!(axC, t, C, label="With Treatment")
    lines!(axC, t, C_no_treatment, linestyle=:dash, label="No Treatment")

    axU = Axis(right[2, 1], title="Control Levels", xticks=50:10:100)
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

    t = result.t .+ 50.0
    sol = result.sol

    names = ["Aβ", "Ca", "τ", "N", "C"]

    fig = Figure(size=(1600, 900))

    for i in 1:5
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1

        ax = Axis(fig[row, col], title=string(names[i], " Levels"), xlabel="Age", xticks=50:10:100)
        lines!(ax, t, sol[i, :])
    end

    return fig
end

end
