"""
    TMJets21a{N, DM<:AbstractDisjointnessMethod} <: AbstractContinuousPost

Set propagation using Taylor models with the function `validated_integ` in `TaylorModels.jl`.
"""
@with_kw struct TMJets21a{N,DM<:AbstractDisjointnessMethod} <: AbstractContinuousPost
    orderQ::Int = DEFAULT_ORDER_Q_TMJETS
    orderT::Int = DEFAULT_ORDER_T_TMJETS
    abstol::N = DEFAULT_ABS_TOL_TMJETS
    maxsteps::Int = DEFAULT_MAX_STEPS_TMJETS
    adaptive::Bool = true
    minabstol::N = Float64(TM._DEF_MINABSTOL)
    absorb::Bool = false
    disjointness::DM = ZonotopeEnclosure()
end

numtype(::TMJets21a{N}) where {N} = N
rsetrep(::TMJets21a{N}) where {N} = TaylorModelReachSet{N}

include("post.jl")
