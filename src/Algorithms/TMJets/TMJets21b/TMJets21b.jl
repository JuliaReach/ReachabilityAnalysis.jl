"""
    TMJets21b{N, DM<:AbstractDisjointnessMethod} <: AbstractContinuousPost

Set propagation using Taylor models with the function `validated_integ` in `TaylorModels.jl`.
"""
@with_kw struct TMJets21b{N,DM<:AbstractDisjointnessMethod} <: AbstractContinuousPost
    orderQ::Int = DEFAULT_ORDER_Q_TMJETS
    orderT::Int = DEFAULT_ORDER_T_TMJETS
    abstol::N = DEFAULT_ABS_TOL_TMJETS
    maxsteps::Int = DEFAULT_MAX_STEPS_TMJETS
    absorb::Bool = false
    adaptive::Bool = true
    minabstol::N = Float64(TM._DEF_MINABSTOL)
    validatesteps::Int = 30
    ε::N = 1e-10
    δ::N = 1e-6
    absorb_steps::Int = 3
    disjointness::DM = ZonotopeEnclosure()
end

numtype(::TMJets21b{N}) where {N} = N
rsetrep(::TMJets21b{N}) where {N} = TaylorModelReachSet{N}

include("post.jl")
