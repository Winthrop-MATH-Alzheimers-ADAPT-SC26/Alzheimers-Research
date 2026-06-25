cd(@__DIR__)
include("modules/Project.jl")
includet("modules/PlottingTools.jl")

using .Project
using .PlottingTools
using LaTeXStrings

# ── Setup ─────────────────────────────────────────────────────────────────────
pvec = build_param_vector(c₃=1.5778463)
tspan = (0.0, 50.0)
umax = (10.0, 1.0e7, 10.0)

u1_treat_weights = calculate_weights(sys, pvec, tspan, (umax[1], 0.0, 0.0))
u2_treat_weights = calculate_weights(sys, pvec, tspan, (0.0, umax[2], 0.0))
u3_treat_weights = calculate_weights(sys, pvec, tspan, (0.0, 0.0, umax[3]))

weight_scales = range(start=1.0, step=0.5, length=4)  # [1.0, 1.5, 2.0, 2.5]

# ── run_scenario ──────────────────────────────────────────────────────────────
function run_scenario(pvec, tspan, umax, u1w, u2w, u3w;
    w3_scale=1.0, w4_scale=1.0, w5_scale=1.0,
    scenario_title="No Title")
    println("Started: $scenario_title")

    p1 = FBSParams(
        u₁max=umax[1], u₂max=0.0, u₃max=0.0,
        w₁=u1w.w1, w₂=u1w.w2, w₃=u1w.w3 * w3_scale, w₄=0.0, w₅=0.0,
        max_iter=10000, relax=0.5, verbose=false
    )
    p2 = FBSParams(
        u₁max=0.0, u₂max=umax[2], u₃max=0.0,
        w₁=u2w.w1, w₂=u2w.w2, w₃=0.0, w₄=u2w.w4 * w4_scale, w₅=0.0,
        max_iter=10000, relax=0.5, verbose=false
    )
    p3 = FBSParams(
        u₁max=0.0, u₂max=0.0, u₃max=umax[3],
        w₁=u3w.w1, w₂=u3w.w2, w₃=0.0, w₄=0.0, w₅=u3w.w5 * w5_scale,
        max_iter=10000, relax=0.5, verbose=false
    )

    r1 = forward_backward_sweep(sys, pvec, tspan, p1, n=52)
    r2 = forward_backward_sweep(sys, pvec, tspan, p2, n=52)
    r3 = forward_backward_sweep(sys, pvec, tspan, p3, n=52)

    println("Finished: $scenario_title")
    return (r1, r2, r3)
end

# ── Baseline ──────────────────────────────────────────────────────────────────
baseline = run_scenario(pvec, tspan, umax,
    u1_treat_weights, u2_treat_weights, u3_treat_weights,
    scenario_title="Baseline")

# ── Vary w3 only ──────────────────────────────────────────────────────────────
scales = weight_scales[2:end]   # [1.5, 2.0, 2.5]

vary_w3 = [run_scenario(pvec, tspan, umax,
    u1_treat_weights, u2_treat_weights, u3_treat_weights,
    w3_scale=s, scenario_title="Varying w₃ s=$s")
           for s in scales]

# ── Collect & plot ────────────────────────────────────────────────────────────
results = (
    baseline=baseline,
    vary_w3=vary_w3,
);

fig = make_plot(WeightSensitivitySingle(), results;
    vary_key=:vary_w3,
    weight_label=L"Vary $w_3$")
