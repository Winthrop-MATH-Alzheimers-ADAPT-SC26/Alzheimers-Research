include("modules/Project.jl")
includet("modules/PlottingTools.jl")
using .Project
using .PlottingTools
using CairoMakie: save

pvec = build_param_vector(c₃=1.5778463)
tspan = (0.0, 50.0)
umax = (10.0, 1.0e7, 10.0)

u1_treat_weights = calculate_weights(sys, pvec, tspan, (umax[1], 0.0, 0.0))
u2_treat_weights = calculate_weights(sys, pvec, tspan, (0.0, umax[2], 0.0))
u3_treat_weights = calculate_weights(sys, pvec, tspan, (0.0, 0.0, umax[3]))

weight_scales = range(start=1.0, step=0.5, length=4)  # [1.0, 1.5, 2.0, 2.5]

function run_scenario(; w3_scale=1.0, w4_scale=1.0, w5_scale=1.0, scenario_title="No Title")
    println("|----- Started Scenario: $scenario_title")
    println("|")

    p1 = FBSParams(
        u₁max=umax[1], u₂max=0.0, u₃max=0.0,
        w₁=u1_treat_weights.w1, w₂=u1_treat_weights.w2, w₃=u1_treat_weights.w3 * w3_scale, w₄=0.0, w₅=0.0,
        max_iter=10000, relax=0.5, verbose=false
    )
    p2 = FBSParams(
        u₁max=0.0, u₂max=umax[2], u₃max=0.0,
        w₁=u2_treat_weights.w1, w₂=u2_treat_weights.w2, w₃=0.0, w₄=u2_treat_weights.w4 * w4_scale, w₅=0.0,
        max_iter=10000, relax=0.5, verbose=false
    )
    p3 = FBSParams(
        u₁max=0.0, u₂max=0.0, u₃max=umax[3],
        w₁=u3_treat_weights.w1, w₂=u3_treat_weights.w2, w₃=0.0, w₄=0.0, w₅=u3_treat_weights.w5 * w5_scale,
        max_iter=10000, relax=0.5, verbose=false
    )

    println("|--- FBS: u₁")
    print("|-- ")
    r1 = forward_backward_sweep(sys, pvec, tspan, p1, n=52)
    println("|--- FBS: u₂")
    print("|-- ")
    r2 = forward_backward_sweep(sys, pvec, tspan, p2, n=52)
    println("|--- FBS: u₃")
    print("|-- ")
    r3 = forward_backward_sweep(sys, pvec, tspan, p3, n=52)
    println("|")
    println("|----- Finished Scenario: $scenario_title")
    return (r1, r2, r3)
end

# Baseline (col 1): scale = 1.0 for all
baseline = run_scenario(scenario_title="Baseline");

# Vary w1 across scales (cols 2-4), w2 and w3 held at baseline
vary_w3 = [run_scenario(w3_scale=s, scenario_title="Varying w₃") for s in weight_scales[2:end]];

# Vary w2 across scales (cols 2-4), w1 and w3 held at baseline
vary_w4 = [run_scenario(w4_scale=s, scenario_title="Varying w₄") for s in weight_scales[2:end]];

# Vary w3 across scales (cols 2-4), w1 and w2 held at baseline
vary_w5 = [run_scenario(w5_scale=s, scenario_title="Varying w₅") for s in weight_scales[2:end]];

results = (
    baseline=baseline,
    vary_w3=vary_w3,
    vary_w4=vary_w4,
    vary_w5=vary_w5,
    scales=weight_scales,
);

fig = make_plot(WeightSensitivity(), results)
#save("figures/weight_sensitivity.png", fig)
