@with_kw struct ASB07{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=CorrectionHullApproximation(order=10, exp_method="base")
    max_order::Int=10
end

include("post.jl")
include("reach.jl")
