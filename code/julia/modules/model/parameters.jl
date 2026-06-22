module Parameters

export build_param_vector

using ..Model: sys

const DEFAULT_PARAMS = (
    Aβ=0.0,
    Ca=100.0,
    τ=0.0,
    N=0.0,
    C=0.0,
    a₁=6.5641,
    a₂=1.5778463,
    b₁=315569260,
    b₂=6311385.2,
    c₁=52.2958,
    c₂=1.78367,
    c₃=0.1,
    d₁=0.07176,
    d₂=0.398406,
    e₁=0.01463,
    e₂=0.008684,
    e₃=0.01992,
    k₁=0.303598,
    k₂=3155692.6,
    k₃=1.8333,
    k₄=0.3588,
    k₅=0.08684,
    σ₁=100.0,
    σ₂=100.0,
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
        sys.σ₁ => p.σ₁,
        sys.σ₂ => p.σ₂,
        sys.R => p.R,
        sys.u₁ => p.u₁,
        sys.u₂ => p.u₂,
        sys.u₃ => p.u₃,
        sys.Aβ => p.Aβ,
        sys.Ca => p.Ca,
        sys.τ => p.τ,
        sys.N => p.N,
        sys.C => p.C
    ]
end

end
