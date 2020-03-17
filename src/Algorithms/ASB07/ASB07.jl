@with_kw struct ASB07{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM # = CorrectionHull(order=10)
    max_order::Int=10
end

include("post.jl")
include("reach.jl")
