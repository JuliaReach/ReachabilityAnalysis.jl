module CorrectionHullMatrixZonotopeModule

using LazySets: ExponentialMap, MatrixZonotope, MatrixZonotopeExp,
                SparsePolynomialZonotope, AbstractZonotope,
                dim, overapproximate
using Reexport: @reexport
using MathematicalSystems: IVP, LinearParametricContinuousSystem,
                           LinearParametricDiscreteSystem,
                           ConstrainedLinearControlParametricContinuousSystem,
                           ConstrainedLinearControlParametricDiscreteSystem,
                           initial_state, state_matrix
using ..DiscretizationModule
using ..ReachabilityAnalysis: IDGenerator, synchronize!, fresh!
using LinearAlgebra: I, norm

export CorrectionHullMatrixZonotope
export overapproximate_continuous_input
export overapproximate_discrete_input
export overapproximate_discrete_input_split

@reexport import ..DiscretizationModule: discretize
@reexport import LazySets: concretize

const LPCS = LinearParametricContinuousSystem
const LPDS = LinearParametricDiscreteSystem
const CLCPCS = ConstrainedLinearControlParametricContinuousSystem
const CLCPDS = ConstrainedLinearControlParametricDiscreteSystem

@inline function _split_emz_overapproximation(MZP,
                                              P::S,
                                              k::Int,
                                              matnorm::Real) where {S<:Union{SparsePolynomialZonotope,
                                                                             AbstractZonotope}}
    tayexp = LazySets.Approximations.taylor_expmap_truncation(MZP, P, k)
    Z = overapproximate(P, Zonotope)
    lagrem = LazySets.Approximations.taylor_expmap_remainder(Z, matnorm, k)
    return tayexp, lagrem
end

function overapproximate_continuous_input(A::MatrixZonotope{N},
                                          B::MatrixZonotope{N},
                                          T::MatrixZonotope{N},
                                          U::SparsePolynomialZonotope{N},
                                          idg::IDGenerator,
                                          taylor_order::Int,
                                          A_norm::N;
                                          Δt::N) where {N}
    fresh!(idg, U)
    BU = overapproximate(B * U, SparsePolynomialZonotope)
    BUT = overapproximate(T * BU, SparsePolynomialZonotope)
    AT = A * T
    AΔt_norm = A_norm * Δt

    poly = LazySets.Approximations.taylor_expmap_truncation(AT, BUT, taylor_order)
    Z = overapproximate(BU, Zonotope)
    rem = scale!(Δt, LazySets.Approximations.taylor_expmap_remainder(Z, AΔt_norm, taylor_order))

    return remove_redundant_generators(minkowski_sum(poly, rem))
end

function overapproximate_discrete_input(A::MatrixZonotope{N},
                                        B::MatrixZonotope{N},
                                        U::SparsePolynomialZonotope{N},
                                        idg::IDGenerator,
                                        taylor_order::Int,
                                        A_norm::N,
                                        t::N) where {N}
    poly, zono = overapproximate_discrete_input_split(A, B, U, idg, taylor_order, A_norm, t)
    return MinkowskiSum(poly, zono)
end

function overapproximate_discrete_input_split(A::MatrixZonotope{N},
                                              B::MatrixZonotope{N},
                                              U::SparsePolynomialZonotope{N},
                                              idg::IDGenerator,
                                              taylor_order::Int,
                                              A_norm::N,
                                              t::N) where {N}

    fresh!(idg, U)
    BUt = scale!(t, overapproximate(B * U, SparsePolynomialZonotope))
    At = scale(t, A)
    At_norm = A_norm * t
    tayexp, lagrem = _split_emz_overapproximation(At, BUt, taylor_order, At_norm)
    return remove_redundant_generators(tayexp), remove_redundant_generators(lagrem)
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
    return CorrectionHullMatrixZonotope{Val{recursive}}(taylor_order, Val(recursive), IDGenerator(0))
end

# Homogeneous case
function discretize(ivp::IVP{<:LPCS,<:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope)
    taylor_order = alg.taylor_order
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)

    N = eltype(X0)

    idg = alg.idg
    synchronize!(idg, A)
    IDₜ = idg(1)
    
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], IDₜ)

    Ω0 = _discretize_CHMZ(A, T, X0, taylor_order, alg.recursive, idg)

    Sdis = LPDS(A)
    return IVP(Sdis, Ω0)
end

# Recursive case
function _discretize_CHMZ(A, T, X0, taylor_order, ::Val{true}, idg)
    _ = idg
    expAT = MatrixZonotopeExp(A * T)
    em = ExponentialMap(expAT, X0)
    return overapproximate(em, SparsePolynomialZonotope, taylor_order)
end

# Non-recursive case
function _discretize_CHMZ(A, T, X0, taylor_order, ::Val{false}, idg::IDGenerator)
    AT = overapproximate(A * T, MatrixZonotope)
    expAT = MatrixZonotopeExp(AT)
    expAT_approx = overapproximate(expAT, MatrixZonotope, taylor_order)
    synchronize!(idg, A, X0)
    fresh!(idg, expAT_approx)
    return overapproximate(expAT_approx * X0, SparsePolynomialZonotope)
end

# Non-homogeneous case
function discretize(ivp::IVP{<:CLCPCS,<:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope)
    taylor_order = alg.taylor_order
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)

    N = eltype(X0)

    idg = alg.idg
    synchronize!(idg, A)
    IDₜ = idg(1)[1]
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], [IDₜ])

    H0 = _discretize_CHMZ(A, T, X0, taylor_order, alg.recursive, idg)

    B = input_matrix(ivp)
    U = inputset(ivp)
    A_norm = norm(A, Inf)
    Pτ0 = overapproximate_continuous_input(A, B, T, U, idg, taylor_order, A_norm; Δt=δ)

    Ω0 = H0 ⊞ Pτ0
    Sdis = CLCPDS(A, B, ivp.s.X, ivp.s.U)
    return InitialValueProblem(Sdis, Ω0)
end

end # module
