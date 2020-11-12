"""
    QINT{N, AM} <: AbstractContinuousPost

Reachability method for one-dimensional quadratic ODEs with non-deterministic inputs.

## Fields

The type fields are:

- `N`   -- number type of the step-size
- `AM`  -- type of the approximation model


## Notes


"""
struct QINT{N, AM} <: AbstractContinuousPost
    Δ::N
    δ::N
    θ::N
    maxiter::Int
    approx_model::AM
end

# convenience constructor using symbols
function QINT(; Δ::N, δ::N, θ::N, maxiter::Int,
               approx_model::AM=Forward(sih=:concrete, exp=:base, setops=:interval)) where {N, AM}
    return QINT(Δ, δ, θ, maxiter, approx_model)
end

step_size(alg::QINT) = alg.δ
numtype(::QINT{N}) where {N} = N
setrep(::QINT{N}) where {N} = Interval{N, IA.Interval{N}}
rsetrep(::QINT{N}) where {N} = ReachSet{N, Interval{N, IA.Interval{N}}}

include("post.jl")
include("reach_homog.jl")
#include("reach_inhomog.jl")
