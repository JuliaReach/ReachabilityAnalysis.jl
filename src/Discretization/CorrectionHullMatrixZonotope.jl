# =======================================================================
# Discretize using the correction hull of the matrix zonotope exponential
# =======================================================================
module CorrectionHullMatrixZonotopeModule

using LazySets
using Reexport

export CorrectionHullMatrixZonotope

@reexport import ..DiscretizationModule: discretize

struct CorrectionHullMatrixZonotope{EM} <: AbstractApproximationModel
    exp::EM
end

# ------------------------------------------------------------
# CorrectionHullMatrixZonotope approximation: Homogeneous case
# ------------------------------------------------------------

function discretize(ivp::IVP{}, δ, alg::CorrectionHullMatrixZonotope)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    n = dim(X0)

    Tₜ = N(0.5) * δ * Matrix(N(1) * I, n, n)
    T = MatrixZonotope(Tₜ, [Tₜ], [100]) #TODO add IDgen
    expAT = MatrixZonotopeExp(A * T)
    
    em = ExponentialMap(expAT, X0)
    Ω0 = overapproximate(em, SparsePolynomialZonotope, taylor_order) #add taylor_order as arg 

    # create result
    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

end # module