"""
    TMJets{N} <: AbstractContinuousPost

Validated integration with Taylor models, based on the algorithm implemented
by Luis Benet and David Sanders in `TalorModels.jl`.

### Fields

- `max_steps` -- (optional, default: `1000`) maximum number of steps in the
                 validated integration ``x' = f(x)``
- `abs_tol`   -- (optional, default: `1e-15`) absolute tolerance
- `orderT`    -- (optional, default: `8`) order of the Taylor model in time
- `orderQ`    -- (optional, default: `2`) order of the Taylor models for jet
                 transport variales
"""
@with_kw struct TMJets{N} <: AbstractContinuousPost
    max_steps::Int=2000
    abs_tol::N=1e-15
    orderT::Int=8
    orderQ::Int=2
end

using TaylorModels: TaylorModelN
using TaylorModels: fp_rpa, remainder

numtype(::TMJets{N}) where {N} = N
rsetrep(::TMJets{N}) where {N} = TaylorModelReachSet{N}

include("post.jl")
include("reach.jl")
include("validated_integ.jl")
