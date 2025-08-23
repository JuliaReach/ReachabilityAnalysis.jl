using LazySets

"""
    HLBS25{N, AM, RM, R} <: AbstractContinuousPost

Implementation of the reachability algorithm for linear systems with parametric
uncertainty using matrix zonotopes by [HuangLBS25](@citet).

### Fields

- `δ`                -- step-size of the discretization
- `approx_model`     -- (optional, default: `CorrectionHullMatrixZonotope`)
                        model for the discretization of the initial value problem
- `max_order`        -- (optional, default: `5`) maximum order of the Taylor
                        series expansion of the matrix exponential
- `taylor_order`     -- (optional, default: `5`) order of the Taylor series
                        expansion of the matrix exponential for each step
- `reduction_method` -- (optional, default: `LazySets.GIR05()`) matrix zonotope
                        order reduction method used
- `recursive`        -- (optional, default: `false`) if `true`, compute the
                        Taylor series expansion of the matrix zonotope
                        exponential map recursively

### Notes

The type parameter `N` is the numeric type, `AM` is the type of the approximation model,
`RM` is the type of the reduction method, and `R` is a `Val{Bool}` indicating if the
recursive method is used.

The `recursive` option is used to compute the Taylor expansion of the matrix
zonotope exponential map. If `recursive == true`, each term of the Taylor expansion
is computed recursively, e.g., ``A^2 P = A (A P)``. If `recursive == false`, the
Taylor expansion is computed by overapproximating the matrix zonotope exponential
map, returning a single matrix that represents the exponential. For more details,
we refer to the documentation of `LazySets.jl`.
"""
struct HLBS25{N,AM,RM,R} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    max_order::Int
    taylor_order::Int
    reduction_method::RM
    recursive::R
end

function HLBS25(; δ::N,
                approx_model::AM=CorrectionHullMatrixZonotope(),
                max_order::Int=5,
                taylor_order::Int=5,
                reduction_method::RM=LazySets.GIR05(),
                recursive::Bool=false) where {N,AM,RM}
    return HLBS25{N,AM,RM,Val{recursive}}(δ, approx_model, max_order, taylor_order,
                                          reduction_method, Val(recursive))
end

step_size(alg::HLBS25) = alg.δ
numtype(::HLBS25{N}) where {N} = N

rsetrep(::HLBS25{N}) where {N} = ReachSet{N,SPZ{N}}

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
