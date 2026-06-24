include("modules/Project.jl")
includet("modules/PlottingTools.jl")
using .Project
using .PlottingTools
using CairoMakie: save

pvec = build_param_vector(c₃=1.5778463)
tspan = (0.0, 50.0)
umax = (10.0, 1.0e7, 10.0)

weights1 = calculate_weights(sys, pvec, tspan, (umax[1], 0.0, 0.0))
weights2 = calculate_weights(sys, pvec, tspan, (0.0, umax[2], 0.0))
weights3 = calculate_weights(sys, pvec, tspan, (0.0, 0.0, umax[3]))

weight_scales = range(start=1.0, step=0.5, length=4)  # [1.0, 1.5, 2.0, 2.5]

function run_scenario(; w1_scale=1.0, w2_scale=1.0, w3_scale=1.0, scenario_title="No Title")
    println("|----- Started Scenario: $scenario_title")
    println("|")

    p1 = FBSParams(
        u₁max=umax[1], u₂max=0.0, u₃max=0.0,
        w₁=weights1.w1, w₂=weights1.w2, w₃=weights1.w3 * w1_scale, w₄=0.0, w₅=0.0,
        max_iter=10000, relax=0.5, verbose=false
    )
    p2 = FBSParams(
        u₁max=0.0, u₂max=umax[2], u₃max=0.0,
        w₁=weights2.w1, w₂=weights2.w2, w₃=0.0, w₄=weights2.w4 * w2_scale, w₅=0.0,
        max_iter=10000, relax=0.5, verbose=false
    )
    p3 = FBSParams(
        u₁max=0.0, u₂max=0.0, u₃max=umax[3],
        w₁=weights3.w1, w₂=weights3.w2, w₃=0.0, w₄=0.0, w₅=weights3.w5 * w3_scale,
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
vary_w1 = [run_scenario(w1_scale=s, scenario_title="Varying w₁") for s in weight_scales[2:end]];

# Vary w2 across scales (cols 2-4), w1 and w3 held at baseline
vary_w2 = [run_scenario(w2_scale=s, scenario_title="Varying w₂") for s in weight_scales[2:end]];

# Vary w3 across scales (cols 2-4), w1 and w2 held at baseline
vary_w3 = [run_scenario(w3_scale=s, scenario_title="Varying w₃") for s in weight_scales[2:end]];

results = (
    baseline=baseline,
    vary_w1=vary_w1,
    vary_w2=vary_w2,
    vary_w3=vary_w3,
    scales=weight_scales,
);

fig = make_plot(WeightSensitivity(), results)
#save("figures/weight_sensitivity.png", fig)
