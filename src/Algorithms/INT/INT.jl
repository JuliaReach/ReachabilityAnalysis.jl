@with_kw struct INT{N, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM=Forward(sih=:concrete, exp=:base, phi2=:base, setops=:Interval)
end

step_size(alg::INT) = alg.δ

include("post.jl")
include("reach.jl")
