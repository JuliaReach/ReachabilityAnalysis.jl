struct QINT{N, AM} <: AbstractContinuousPost
    Δ::N
    δ::N
    θ::N
    maxiter::Int
    approx_model::AM
end

# convenience constructor using symbols
function QINT(; Δ::N, δ::N, θ::N, maxiter::Int
               approx_model::AM=Forward(sih=:concrete, exp=:base, setops=:interval)) where {N, AM}
    return INT(Δ, δ, θ, maxiter, approx_model)
end

step_size(alg::INT) = alg.δ
numtype(::INT{N}) where {N} = N
setrep(::INT{N}) where {N} = Interval{N, IA.Interval{N}}
rsetrep(::INT{N}) where {N} = ReachSet{N, Interval{N, IA.Interval{N}}}

include("post.jl")
include("reach_homog.jl")
#include("reach_inhomog.jl")
