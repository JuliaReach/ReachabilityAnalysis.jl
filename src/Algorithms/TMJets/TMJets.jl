"""
    TMJets{N, DM<:AbstractDisjointnessMethod} <: AbstractContinuousPost

Validated integration using Taylor models.

### Fields

- `max_steps`    -- (optional, default: `2000`) maximum number of steps in the
                    validated integration ``x' = f(x)``
- `abs_tol`      -- (optional, default: `1e-10`) absolute tolerance
- `orderT`       -- (optional, default: `8`) order of the Taylor model in time
- `orderQ`       -- (optional, default: `2`) order of the Taylor models for jet transport variables
- `disjointness` -- (optional, default: `ZonotopeEnclosure()`) defines the method to
                    perform the disjointness check between the taylor model flowpipe and the invariant
- `adaptive`     -- (optional, default: `true`) if `true`, try decreasing the absolute
                    tolerance each time step validation fails, until `min_abs_tol` is reached
- `min_abs_tol`  -- (optional, default: `1e-29`) minimum absolute tolerance for the adaptive algorithm

### Notes

The argument `disjointness` allows to control how are disjointness checks
computed, in the case where the invariant is not universal. In particular,
`ZonotopeEnclosure()` pre-processes the taylor model with a zonotopic overapproximation,
then performs the disjointness check with that zonotope and the invariant.
For other options, see the documentation of `AbstractDisjointnessMethod`.

This algorithm is an adaptation of the implementation in `TaylorModels.jl`
(see copyright license in `src/Algorithms/TMJets/reach.jl`). The package
`TaylorIntegration.jl` is used for jet-transport of ODEs using the Taylor method,
and `TaylorSeries.jl` is used to work with truncated Taylor series.
"""
@with_kw struct TMJets{N, DM<:AbstractDisjointnessMethod} <: AbstractContinuousPost
    max_steps::Int=DEFAULT_MAX_STEPS_TMJETS
    abs_tol::N=DEFAULT_ABS_TOL_TMJETS
    orderT::Int=DEFAULT_ORDER_T_TMJETS
    orderQ::Int=DEFAULT_ORDER_Q_TMJETS
    disjointness::DM=ZonotopeEnclosure()
    adaptive::Bool=true
    min_abs_tol::N=1e-29
end

const DEFAULT_MAX_STEPS_TMJETS = 2000
const DEFAULT_ABS_TOL_TMJETS = 1e-10
const DEFAULT_ORDER_T_TMJETS = 8
const DEFAULT_ORDER_Q_TMJETS = 2

using TaylorModels: TaylorModelN
using TaylorModels: fp_rpa, remainder

numtype(::TMJets{N}) where {N} = N
rsetrep(::TMJets{N}) where {N} = TaylorModelReachSet{N}

# TODO refactor
_is_intersection_empty(Ri::TaylorModelReachSet, X::Universe, disjointness::AbstractDisjointnessMethod) = false
_is_intersection_empty(Ri::TaylorModelReachSet, X::Universe, disjointness::ZonotopeEnclosure) = false

include("post.jl")
include("reach.jl")
