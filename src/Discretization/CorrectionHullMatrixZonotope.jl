module CorrectionHullMatrixZonotopeModule

using LazySets: ExponentialMap, MatrixZonotope, MatrixZonotopeExp,
                SparsePolynomialZonotope, dim, indexvector, ngens,
                overapproximate
using Reexport: @reexport
using MathematicalSystems: IVP, LinearParametricContinuousSystem,
                           LinearParametricDiscreteSystem, initial_state,
                           state_matrix
using ..DiscretizationModule
using LinearAlgebra: I

export CorrectionHullMatrixZonotope

@reexport import ..DiscretizationModule: discretize

const LPCS = LinearParametricContinuousSystem
const LPDS = LinearParametricDiscreteSystem

"""
    CorrectionHullMatrixZonotope{R} <: AbstractApproximationModel

Implementation of the discretization algorithm for linear systems with parametric
uncertainty using matrix zonotopes by [HuangLBS25](@citet).

### Fields

- `taylor_order`     -- (optional, default: `5`) order of the Taylor series
                        expansion of the matrix exponential for each step
- `recursive`        -- (optional, default: `false`) if `true`, compute the
                        Taylor series expansion of the matrix zonotope
                        exponential map recursively

### Notes

The `recursive` option is used to compute the Taylor expansion of the matrix zonotope exponential map.
If `recursive == true`, each term of the Taylor expansion is computed recursively (e.g., ``A^2 P = A (A P)``).

If `recursive == false`, the Taylor expansion is computed by overapproximating the matrix zonotope exponential
map, producing a single matrix that represents the exponential.

"""
struct CorrectionHullMatrixZonotope{R} <: AbstractApproximationModel
    taylor_order::Int
    recursive::R
end

function CorrectionHullMatrixZonotope(; taylor_order::Int=5, recursive::Bool=false)
    return CorrectionHullMatrixZonotope{Val{recursive}}(taylor_order, Val(recursive))
end

function discretize(ivp::IVP{<:LPCS,<:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope)
    taylor_order = alg.taylor_order
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)

    N = eltype(X0)

    IDₜ = ngens(A) > 0 ? maximum(indexvector(A)) + 1 : 1
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], [IDₜ])

    Ω0 = _discretize_CHMZ(A, T, X0, taylor_order, alg.recursive)

    Sdis = LPDS(A)
    return IVP(Sdis, Ω0)
end

# Recursive case
function _discretize_CHMZ(A, T, X0, taylor_order, ::Val{true})
    expAT = MatrixZonotopeExp(A * T)
    em = ExponentialMap(expAT, X0)
    return overapproximate(em, SparsePolynomialZonotope, taylor_order)
end

# Non-recursive case
function _discretize_CHMZ(A, T, X0, taylor_order, ::Val{false})
    AT = overapproximate(A * T, MatrixZonotope)
    expAT = MatrixZonotopeExp(AT)
    expAT_approx = overapproximate(expAT, MatrixZonotope, taylor_order)
    return overapproximate(expAT_approx * X0, SparsePolynomialZonotope)
end

end # module
