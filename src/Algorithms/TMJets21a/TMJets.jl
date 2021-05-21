const TMJets = TMJets21b

@with_kw struct TMJets20a{N, DM<:AbstractDisjointnessMethod} <: AbstractContinuousPost
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
using TaylorModels: fp_rpa, remainder, initialize!

numtype(::TMJets20{N}) where {N} = N
rsetrep(::TMJets20{N}) where {N} = TaylorModelReachSet{N}

include("post.jl")
include("reach.jl")
