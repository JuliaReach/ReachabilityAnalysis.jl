module CorrectionHullMatrixZonotopeModule

using LazySets, Reexport
using MathematicalSystems
using ..DiscretizationModule
using LinearAlgebra: I

export CorrectionHullMatrixZonotope

@reexport import ..DiscretizationModule: discretize

const LPCS = LinearParametricContinuousSystem
const LPDS = LinearParametricDiscreteSystem

struct CorrectionHullMatrixZonotope{N, R} <: AbstractApproximationModel
    taylor_order::Int
    recursive::R
end

function CorrectionHullMatrixZonotope(N::Type=Float64; taylor_order::Int=5, recursive::Bool=true)
    return CorrectionHullMatrixZonotope{N, Val{recursive}}(taylor_order, Val(recursive))
end

function discretize(ivp::IVP{<:LPCS, <:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope{N, Val{true}}) where {N}
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

    X = stateset(ivp)
    Sdis = LPDS(A)
    return InitialValueProblem(Sdis, Ω0)
end

function discretize(ivp::IVP{<:LPCS, <:SparsePolynomialZonotope}, δ,
                    alg::CorrectionHullMatrixZonotope{N, Val{false}}) where {N}
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

    X = stateset(ivp)
    Sdis = LPDS(A)
    return InitialValueProblem(Sdis, Ω0)
end

end # module
