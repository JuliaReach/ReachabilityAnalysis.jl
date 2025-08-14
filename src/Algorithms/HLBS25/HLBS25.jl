using LazySets

struct HLBS25{N,AM,RM,R} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    max_order::Int
    taylor_order::Int
    reduction_method::RM
    recursive::R
end

function HLBS25(; δ::N,
                approx_model::AM=CorrectionHullMatrixZonotope(; order=5),
                max_order::Int=5,
                taylor_order::Int=5,
                reduction_method::RM=LazySets.GIR05(),
                recursive::Bool=false) where {N,AM,RM}
    return HLBS25{N,AM,RM,Val{recursive}}(δ, approx_model, max_order, taylor_order,
                                          reduction_method, Val(recursive))
end

HLBS25(δ::N; kwargs...) where {N} = HLBS25(; δ=δ, kwargs...)

step_size(alg::HLBS25) = alg.δ
numtype(::HLBS25{N}) where {N} = N

# Reach-set representation; specialize on the recursive flag if needed
function rsetrep(::HLBS25{N,AM,RM,Val{false}}) where {N,AM,RM}
    return ReachSet{N,
                    SparsePolynomialZonotope{N,Vector{N},Matrix{N},Matrix{N},Matrix{Int},
                                             Vector{Int}}}
end

function rsetrep(::HLBS25{N,AM,RM,Val{true}}) where {N,AM,RM}
    return ReachSet{N,
                    SparsePolynomialZonotope{N,Vector{N},Matrix{N},Matrix{N},Matrix{Int},
                                             Vector{Int}}}
end

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
