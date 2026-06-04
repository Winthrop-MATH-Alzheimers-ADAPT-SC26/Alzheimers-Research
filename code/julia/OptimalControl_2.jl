using Plots
using Plots.PlotMeasures
using PlotThemes
include("FBS.jl")
include("WeightCalculations.jl")

# Plot.jl Defaults
theme(
    :wong,                            # Colorblind-friendly, publication-ready colors
    framestyle=:box,                # Boxed axes instead of floating T-axes
    grid=true,                      # Faint background grid lines for readability
    gridalpha=0.3,                  # Make grid lines subtle
    linewidth=2,                    # Thicker, more visible lines
    markersize=6,                   # Properly scaled scatter markers
    margin=5mm,                     # Ensure axis labels aren't cut off
    dpi=100                         # Crisp resolution for exports
)

# Set u₁, u₂, u₃ maxes
const u1max = 0.12012
const u2max = 0.12012
const u3max = 0.12012
# Set Timespan
const ti = 0.0
const te = 100.0

# Unpack
tspan = (ti, te)
umax = (u1max, u2max, u3max)

# Calculate Weights
weights = calculate_weights(umax=umax, tspan=tspan)

result = forward_backward_sweep(
    x0=[0.0, 100.0, 0.0, 0.0, 0.0],
    tspan=tspan,
    umax=umax,
    weights=weights,
    max_iter=500,
    tol=1e-6
)

# Unpack
ts = result.t;
u1, u2, u3 = result.controls.u1, result.controls.u2, result.controls.u3;
fwd_sol = result.states;

# Plot states
p1 = plot(fwd_sol, idxs=1, xlabel="years", label="Aβ");
p2 = plot!(fwd_sol, idxs=2, xlabel="years", label="Ca");
p3 = plot!(fwd_sol, idxs=3, xlabel="years", label="τ");
p4 = plot!(fwd_sol, idxs=4, title="Aβ, Ca, τ, N States", xlabel="years", label="N");
p5 = plot(fwd_sol, idxs=5, title="C State", xlabel="years", label="C");
# Plot controls
p6 = plot(ts, u1, title="u₁, u₂, u₃ Treatments", xlabel="years", label="u₁");
p7 = plot!(ts, u2, label="u₂");
p8 = plot!(ts, u3, label="u₃");
display(plot(p1, p5, p6, layout=(3, 1), size=(1200, 900)))

println("Final C  = ", fwd_sol[5, end])
println("Final N  = ", fwd_sol[4, end])
