module PlottingTools

export ControlsPlot,
    ControlsPlotSeperate,
    ControlsPlotSeperate2,
    ControlsPlotSeperate3,
    SolPlotCombined,
    SolPlotHorizontal,
    make_plot

using CairoMakie

abstract type PlotType end

struct ControlsPlot <: PlotType end
struct ControlsPlotSeperate <: PlotType end
struct ControlsPlotSeperate2 <: PlotType end
struct ControlsPlotSeperate3 <: PlotType end
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
    t = r1.t
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

    axc1 = Axis(left[1, 1], title=L"\text{C with only } u_1 \text{ Treatment}", xticks=0:10:50)
    lines!(axc1, t, C1, label="With Treatment")
    lines!(axc1, t, C1_base, linestyle=:dash, label="No Treatment")
    axislegend(axc1)

    axc2 = Axis(left[2, 1], title=L"\text{C with only } u_2 \text{ Treatment}", xticks=0:10:50)
    lines!(axc2, t, C2, label="With Treatment")
    lines!(axc2, t, C2_base, linestyle=:dash, label="No Treatment")
    axislegend(axc2)

    axc3 = Axis(left[3, 1], title=L"\text{C with only } u_3 \text{ Treatment}", xticks=0:10:50)
    lines!(axc3, t, C3, label="With Treatment")
    lines!(axc3, t, C3_base, linestyle=:dash, label="No Treatment")
    axislegend(axc3)

    axc4l = Axis(right[1, 1], title="Each Treatment", xticks=0:10:50, ylabel=L"u_1 \text{ & } u_2")
    axc4r = Axis(right[1, 1], xticks=0:10:50, yaxisposition=:right, ylabel=L"u_2")
    l1 = lines!(axc4l, t, u1, label=L"u_1")
    l3 = lines!(axc4l, t, u3, label=L"u_3")
    l2 = lines!(axc4r, t, u2, label=L"u_2")

    Legend(
        right[1, 1],
        [l1, l2, l3],
        [L"u_1", L"u_2", L"u_3"],
        tellwidth=false,
        tellheight=false,
        halign=:right,
        valign=:top
    )
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)

    return fig
end

function make_plot(::ControlsPlotSeperate2, r1, r2, r3)
    t = r1.t
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

    axc1 = Axis(left[1, 1], title=L"\text{C with only } u_1 \text{ Treatment}", xticks=0:10:50)
    lines!(axc1, t, C1, label="With Treatment")
    lines!(axc1, t, C1_base, linestyle=:dash, label="No Treatment")
    axislegend(axc1)

    axc2 = Axis(left[2, 1], title=L"\text{C with only } u_2 \text{ Treatment}", xticks=0:10:50)
    lines!(axc2, t, C2, label="With Treatment")
    lines!(axc2, t, C2_base, linestyle=:dash, label="No Treatment")
    axislegend(axc2)

    axc3 = Axis(left[3, 1], title=L"\text{C with only } u_3 \text{ Treatment}", xticks=0:10:50)
    lines!(axc3, t, C3, label="With Treatment")
    lines!(axc3, t, C3_base, linestyle=:dash, label="No Treatment")
    axislegend(axc3)

    axc4 = Axis(right[1, 1], title=L"u_1 \text{ Treatment}", xticks=0:10:50)
    lines!(axc4, t, u1, label=L"u_1")
    axislegend(axc4)

    axc5 = Axis(right[2, 1], title=L"u_2 \text{ Treatment}", xticks=0:10:50)
    lines!(axc5, t, u2, label=L"u_2")
    axislegend(axc5)

    axc6 = Axis(right[3, 1], title=L"u_3 \text{ Treatment}", xticks=0:10:50)
    lines!(axc6, t, u3, label=L"u_3")
    axislegend(axc6)

    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)
    return fig
end


function make_plot(::ControlsPlotSeperate3, r1, r2, r3)
    t = r1.t
    sol1, sol2, sol3 = r1.sol, r2.sol, r3.sol

    N_base = r1.baseline_sol[4, :]
    N_u1 = sol1[4, :]
    N_u2 = sol2[4, :]
    N_u3 = sol3[4, :]

    C_base = r1.baseline_sol[5, :]
    C_u1 = sol1[5, :]
    C_u2 = sol2[5, :]
    C_u3 = sol3[5, :]

    u1 = r1.controls.u1
    u2 = r2.controls.u2
    u3 = r3.controls.u3

    fig = Figure(size=(1600, 900))
    left = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()

    # --- Left: 2x1 grid of state comparisons --- #

    axn = Axis(left[1, 1], title=L"\text{Neuron Loss } (N) \text{ - Control Scenario}", xticks=0:10:50, xlabel="Years After Age 50")
    lines!(axn, t, N_base, linestyle=:dash, label="No Treatment")
    lines!(axn, t, N_u1, label=L"N \text{ with only } u_1")
    lines!(axn, t, N_u2, label=L"N \text{ with only } u_2")
    lines!(axn, t, N_u3, label=L"N \text{ with only } u_3")
    axislegend(axn, position=:lt)

    axc = Axis(left[2, 1], title=L"\text{Cognitive Decline } (C) \text{ - Control Scenario}", xticks=0:10:50, xlabel="Years After Age 50")
    lines!(axc, t, C_base, linestyle=:dash, label="No Treatment")
    lines!(axc, t, C_u1, label=L"C \text{ with only } u_1")
    lines!(axc, t, C_u2, label=L"C \text{ with only } u_2")
    lines!(axc, t, C_u3, label=L"C \text{ with only } u_3")
    axislegend(axc, position=:lt)

    # --- Right: 3x1 grid of control profiles over time --- #

    axu1 = Axis(right[1, 1], title=L"u_1 \text{ Treatment}", xticks=0:10:50, xlabel="Years After Age 50")
    lines!(axu1, t, u1, label=L"u_1", color=Cycled(2))
    axislegend(axu1)

    axu2 = Axis(right[2, 1], title=L"u_2 \text{ Treatment}", xticks=0:10:50, xlabel="Years After Age 50")
    lines!(axu2, t, u2, label=L"u_2", color=Cycled(3))
    axislegend(axu2)

    axu3 = Axis(right[3, 1], title=L"u_3 \text{ Treatment}", xticks=0:10:50, xlabel="Years After Age 50")
    lines!(axu3, t, u3, label=L"u_3", color=Cycled(4))
    axislegend(axu3)

    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)
    return fig
end

# Plot controls alongside solution
function make_plot(::ControlsPlot, result)

    t = result.t
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

    Axis(left[1, 1], title=L"Aβ\text{ Levels}", xticks=0:10:50)
    lines!(t, Aβ)
    Axis(left[1, 2], title=L"Ca\text{ Levels}", xticks=0:10:50)
    lines!(t, Ca)
    Axis(left[2, 1], title=L"τ\text{ Levels}", xticks=0:10:50)
    lines!(t, τ)
    Axis(left[2, 2], title=L"N\text{ Levels}", xticks=0:10:50)
    lines!(t, N)

    axC = Axis(right[1, 1], title=L"C\text{ Levels}", xticks=0:10:50)
    lines!(axC, t, C, label="With Treatment")
    lines!(axC, t, C_no_treatment, linestyle=:dash, label="No Treatment")

    axc4l = Axis(right[2, 1], title="Each Treatment", xticks=0:10:50, ylabel=L"u_1 \text{ & } u_2")
    axc4r = Axis(right[2, 1], xticks=0:10:50, yaxisposition=:right, ylabel=L"u_2")
    l1 = lines!(axc4l, t, u1, label=L"u_1")
    l3 = lines!(axc4l, t, u3, label=L"u_3")
    l2 = lines!(axc4r, t, u2, label=L"u_2")

    Legend(
        right[1, 1],
        [l1, l2, l3],
        [L"u_1", L"u_2", L"u_3"],
        tellwidth=false,
        tellheight=false,
        halign=:right,
        valign=:top
    )

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

    lines!(ax, t, sol[1, :], label=L"Aβ")
    lines!(ax, t, sol[2, :], label=L"Ca")
    lines!(ax, t, sol[3, :], label=L"τ")
    lines!(ax, t, sol[4, :], label=L"N")
    lines!(ax, t, sol[5, :], label=L"C")

    axislegend(ax)
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)

    return fig
end

# Individual state plots
function make_plot(::SolPlotHorizontal, result)

    t = result.t
    sol = result.sol

    names = [L"Aβ", L"Ca", L"τ", L"N", L"C"]

    fig = Figure(size=(1600, 900))

    for i in 1:5
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1

        ax = Axis(fig[row, col], title=string(names[i], " Levels"), xlabel="Years After Age 50", xticks=0:10:50)
        lines!(ax, t, sol[i, :])
    end
    colgap!(fig.layout, 15)
    rowgap!(fig.layout, 15)

    return fig
end

end
