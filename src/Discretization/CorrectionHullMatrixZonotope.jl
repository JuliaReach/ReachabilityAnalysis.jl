module CorrectionHullMatrixZonotopeModule

using LazySets, Reexport
using MathematicalSystems
using ..DiscretizationModule
using LinearAlgebra: I

export CorrectionHullMatrixZonotope
export ExactSum

@reexport import ..DiscretizationModule: discretize
@reexport import LazySets: concretize

const LPCS = LinearParametricContinuousSystem
const LPDS = LinearParametricDiscreteSystem

#TODO: move to LazySets
struct ExactSum{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function ExactSum(X::LazySet{N}, Y::LazySet{N}) where {N}
        @assert dim(X) == dim(Y) "sets in a exact sum must have the same dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

function concretize(ES::ExactSum)
    exact_sum(ES.X, ES.Y)
end

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
    idg::IDGenerator
end

function CorrectionHullMatrixZonotope(; taylor_order::Int=5, recursive::Bool=false)
    return CorrectionHullMatrixZonotope{Val{recursive}}(taylor_order, Val(recursive))
end

# Homogeneous case
function discretize(ivp::IVP{<:LPCS,<:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope)
    taylor_order = alg.taylor_order
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)

    N = eltype(X0)

    synchronize!(idg, A)
    IDₜ = idg(1) #IDₜ = ngens(A) > 0 ? maximum(indexvector(A)) + 1 : 1
    
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], IDₜ)

    Ω0 = _discretize_CHMZ(A, T, X0, taylor_order, alg.recursive)

    Sdis = LPDS(A)
    return InitialValueProblem(Sdis, Ω0)
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

# Non-homogeneous case
function discretize(ivp::IVP{<:LCPCS,<:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope)
    taylor_order = alg.taylor_order
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)

    N = eltype(X0)

    IDₜ = ngens(A) > 0 ? maximum(indexvector(A)) + 1 : 1
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], [IDₜ])

    H0 = _discretize_CHMZ(A, T, X0, taylor_order, alg.recursive)

    B = input_matrix(ivp)
    Pτ0 = overapproximate_continuous_input() #TODO implemement 

    Ω0 = ExactSum(H0, Pτ0)
    Sdis = LCPDS(A, B)
    return InitialValueProblem(Sdis, Ω0)
end

end # module
