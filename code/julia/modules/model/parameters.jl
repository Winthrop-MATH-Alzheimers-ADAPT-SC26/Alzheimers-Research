module Parameters

export build_param_vector

using ..Model: sys

const DEFAULT_PARAMS = (
    a₁=65641.0,
    a₂=15778.463,
    b₁=315569260.0,
    b₂=6311385.2,
    c₁=52.2958,
    c₂=7.8,
    c₃=0.0,
    d₁=0.07176,
    d₂=0.3588,
    e₁=0.0,
    e₂=86.84,
    e₃=199.16,
    k₁=3035.98,
    k₂=3155692.6,
    k₃=10.9,
    k₄=0.003588,
    k₅=0.0,
    σ=0.00027,
    R=1.0,
    u₁=0.0,
    u₂=0.0,
    u₃=0.0,
)

function build_param_vector(; kwargs...)
    p = merge(DEFAULT_PARAMS, (; kwargs...))

    return [
        sys.a₁ => p.a₁,
        sys.a₂ => p.a₂,
        sys.b₁ => p.b₁,
        sys.b₂ => p.b₂,
        sys.c₁ => p.c₁,
        sys.c₂ => p.c₂,
        sys.c₃ => p.c₃,
        sys.d₁ => p.d₁,
        sys.d₂ => p.d₂,
        sys.e₁ => p.e₁,
        sys.e₂ => p.e₂,
        sys.e₃ => p.e₃,
        sys.k₁ => p.k₁,
        sys.k₂ => p.k₂,
        sys.k₃ => p.k₃,
        sys.k₄ => p.k₄,
        sys.k₅ => p.k₅,
        sys.σ => p.σ,
        sys.R => p.R,
        sys.u₁ => p.u₁,
        sys.u₂ => p.u₂,
        sys.u₃ => p.u₃,
    ]
end

end
