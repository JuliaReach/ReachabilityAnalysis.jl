@with_kw struct ASB07{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=CorrectionHull(order=10, exp=:base)
    max_order::Int=10
    static::Bool=false
    recursive::Bool=false
end

step_size(alg::ASB07) = alg.δ
numtype(::ASB07{N}) where {N} = N
function rsetrep(::ASB07{N}) where {N}
    if static
        error("not implemented")
    else
        RT = ReachSet{N, Zonotope{N, Vector{N}, Matrix{N}}}
    end
    return RT
end

include("post.jl")
include("reach.jl")
