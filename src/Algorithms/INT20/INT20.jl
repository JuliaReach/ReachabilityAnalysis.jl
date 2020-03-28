@with_kw struct INT20{N, AM} <: AbstractContinuousPost
    δ::N
    # TODO: use set-operations option
    approx_model::AM=Forward(sih_method=:concrete, exp_method=:base,
                             phi2_method=:base, set_operations=:Interval)
end

step_size(alg::INT20) = alg.δ

include("post.jl")
include("reach.jl")
