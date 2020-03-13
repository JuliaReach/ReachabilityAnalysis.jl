"""
    TMJets{N} <: AbstractContinuousPost

Validated integration with Taylor models, based on the algorithm implemented
by Luis Benet and David Sanders in `TalorModels.jl`.

### Fields

- `max_steps` -- maximum number of steps in the validated integration
- `abs_tol`   -- absolute tolerance
- `orderT`    -- order of the Taylor model in time
- `orderQ`    -- order of the Taylor models for jet transport variales
"""
@with_kw struct TMJets{N} <: AbstractContinuousPost
    max_steps::Int=1000
    abs_tol::N=1e-10
    orderT::Int=8
    orderQ::Int=2
    # setrep::ST=Zonotope{Float64, Vector{Float64}, Matrix{Float64}}
    # add setrep as a type?
end

using TaylorModels: TaylorModelN
using TaylorModels: fp_rpa, remainder

include("post.jl")
include("reach.jl")
include("validated_integ.jl")
