module CorrectionHullMatrixZonotopeModule

using LazySets, Reexport
using MathematicalSystems
using ..DiscretizationModule
using LinearAlgebra: I

export CorrectionHullMatrixZonotope

@reexport import ..DiscretizationModule: discretize

const LPCS = LinearParametricContinuousSystem
const LPDS = LinearParametricDiscreteSystem

"""
    CorrectionHullMatrixZonotope{N, R} <: AbstractApproximationModel

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
struct CorrectionHullMatrixZonotope{N,R} <: AbstractApproximationModel
    taylor_order::Int
    recursive::R
    tol::Real
    norm::Real
end

function CorrectionHullMatrixZonotope(N::Type=Float64; taylor_order::Int=5, recursive::Bool=true)
    return CorrectionHullMatrixZonotope{N,Val{recursive}}(taylor_order, Val(recursive))
end

function discretize(ivp::IVP{<:LPCS,<:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope{N,Val{true}}) where {N}
    taylor_order = alg.taylor_order
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)

    IDₜ = ngens(A) > 0 ? maximum(indexvector(A)) + 1 : 1
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], [IDₜ])
    expAT = MatrixZonotopeExp(A * T)

    # recursive logic
    em = ExponentialMap(expAT, X0)
    Ω0 = overapproximate(em, SparsePolynomialZonotope, taylor_order)

    Sdis = LPDS(A)
    return InitialValueProblem(Sdis, Ω0)
end

function discretize(ivp::IVP{<:LPCS,<:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope{N,Val{false}}) where {N}
    taylor_order = alg.taylor_order
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)
    n = dim(X0)

    IDₜ = ngens(A) > 0 ? maximum(indexvector(A)) + 1 : 1
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], [IDₜ])
    AT = overapproximate(A * T, MatrixZonotope)
    expAT = MatrixZonotopeExp(AT)

    # non-recursive logic
    expAT_approx = overapproximate(expAT, MatrixZonotope, taylor_order)
    Ω0 = overapproximate(expAT_approx * X0, SparsePolynomialZonotope)

    Sdis = LPDS(A)
    return InitialValueProblem(Sdis, Ω0)
end

end # module
