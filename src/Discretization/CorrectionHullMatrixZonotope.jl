module CorrectionHullMatrixZonotopeModule

using LazySets
using Reexport

export CorrectionHullMatrixZonotope

@reexport import ..DiscretizationModule: discretize

struct CorrectionHullMatrixZonotope{RM,R} <: AbstractApproximationModel
    taylor_order::Int
    max_order::Int
    reduction_method::RM
    recursive::R
end

function CorrectionHullMatrixZonotope(; taylor_order::Int=5,
                                      max_order::Int=5,
                                      recursive::Bool=false)
    return CorrectionHullMatrixZonotope{Val{recursive}}(taylor_order, max_order, Val(recursive))
end

function discretize(ivp::IVP{<:AbstractContinuousSystem}, δ,
                    alg::CorrectionHullMatrixZonotope{Val{true}})
    @unpack taylor_order, max_order = alg
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
    Sdis = LinearParametricDiscreteSystem(A, X)
    return InitialValueProblem(Sdis, Ω0)
end

function discretize(ivp::IVP{<:AbstractContinuousSystem}, δ,
                    alg::CorrectionHullMatrixZonotope{Val{false}})
    @unpack taylor_order, max_order = alg
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)
    n = dim(X0)

    IDₜ = ngens(A) > 0 ? maximum(indexvector(A)) + 1 : 1
    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], [IDₜ])
    expAT = MatrixZonotopeExp(A * T)

    # non-recursive logic
    expAT_approx = overapproximate(expAT, MatrixZonotope, taylor_order)
    Ω0 = overapproximate(expAT_approx * X0, SparsePolynomialZonotope, taylor_order)

    X = stateset(ivp)
    Sdis = LinearParametricDiscreteSystem(A, X)
    return InitialValueProblem(Sdis, Ω0)
end

end # module
