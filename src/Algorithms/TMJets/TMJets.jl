"""
    TMJets{N, DM<:AbstractIntersectionMethod} <: AbstractContinuousPost

Validated integration with Taylor models, based on the algorithm implemented
by Luis Benet and David Sanders in `TalorModels.jl`.

### Fields

- `max_steps` -- (optional, default: `2000`) maximum number of steps in the
                 validated integration ``x' = f(x)``
- `abs_tol`   -- (optional, default: `1e-10`) absolute tolerance
- `orderT`    -- (optional, default: `8`) order of the Taylor model in time
- `orderQ`    -- (optional, default: `2`) order of the Taylor models for jet transport variables
- `intersection_method` -- (optional, default: `ZonotopeEnclosure()`) defines the method to
                           compute the intersection of the taylor model flowpipe with the invariant
- `adaptive`  -- (optional, default: `true`) if `true`, try decreasing the absolute
                 tolerance each time step validation fails, until `min_abs_tol` is reached
- `min_abs_tol` -- (optional, default: `1e-29`) minimum absolute tolerance for the adaptive algorithm

### Notes

The argument `intersection_method` allows to control how are disjointness checks
computed, in the case where the invariant is not universal. In particular,
`ZonotopeEnclosure()` computes an overapproximation of the taylor model using a
zonotope, then performs the disjointness check with that zonotope and the invariant.
For other methods, see the documentation of `AbstractIntersectionMethod`.
"""
@with_kw struct TMJets{N, DM<:AbstractIntersectionMethod} <: AbstractContinuousPost
    max_steps::Int=2000
    abs_tol::N=1e-10
    orderT::Int=8
    orderQ::Int=2
    intersection_method::AbstractIntersectionMethod=ZonotopeEnclosure()
    adaptive::Bool=true
    min_abs_tol::N=1e-29
end

using TaylorModels: TaylorModelN
using TaylorModels: fp_rpa, remainder

numtype(::TMJets{N}) where {N} = N
rsetrep(::TMJets{N}) where {N} = TaylorModelReachSet{N}

include("post.jl")
include("reach.jl")
include("validated_integ.jl")
