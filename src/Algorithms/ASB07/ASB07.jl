@with_kw struct ASB07{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=CorrectionHull(order=10, exp=:base)
    max_order::Int=10
end

step_size(alg::ASB07) = alg.δ
numtype(::ASB07{N}) where {N} = N

include("post.jl")
include("reach.jl")
