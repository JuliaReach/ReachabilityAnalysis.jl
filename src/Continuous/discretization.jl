abstract type AbstractApproximationModel end
const AAModel = AbstractApproximationModel

export ForwardApproximation,
       BackwardApproximation,
       DiscreteApproximation,
       CorrectionHullApproximation

@with_kw struct ForwardApproximation <: AbstractApproximationModel
    exp_method::String="base"
    set_operations::String="lazy"
    sih_method::String="concrete"
    phi2_method::String="base"
end

@with_kw struct BackwardApproximation <: AbstractApproximationModel
    exp_method::String="base"
    set_operations::String="lazy"
    sih_method::String="concrete"
    phi2_method::String="base"
end

# no bloating
struct DiscreteApproximation <: AbstractApproximationModel
#
end

@with_kw struct CorrectionHullApproximation <: AbstractApproximationModel
   order::Int=10
   exp_method::String="base"
end

function _default_approximation_model(ivp::IVP{<:AbstractContinuousSystem})
    return ForwardApproximation()
end

# homogeneous case
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ::Float64, alg::ForwardApproximation)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    ϕ = _exp(A, δ, alg.exp_method)
    A_abs = _elementwise_abs(A)
    Phi2A_abs = Φ₂(A_abs, δ, alg.phi2_method)

    # "forward" algorithm, uses E⁺
    @assert alg.sih_method == "concrete"
    # TODO : specialize, add option to compute the concrete linear map
    Einit = symmetric_interval_hull(Phi2A_abs * symmetric_interval_hull((A * A) * X0))

    Ω0 = ConvexHull(X0, ϕ * X0 ⊕ Einit)
    X = stateset(ivp)
    Sdiscr = ConstrainedLinearDiscreteSystem(ϕ, X)
    return InitialValueProblem(Sdiscr, Ω0)
end

# inhomogeneous case
function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ::Float64, alg::ForwardApproximation)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)
    U = next_set(inputset(ivp), 1)
    ϕ = _exp(A, δ, alg.exp_method)
    A_abs = _elementwise_abs(A)
    Phi2A_abs = Φ₂(A_abs, δ, alg.phi2_method)

    @assert alg.sih_method == "concrete"
    # TODO : specialize, add option to compute the concrete linear map
    Einit = symmetric_interval_hull(Phi2A_abs * symmetric_interval_hull((A * A) * X0))

    Eψ0 = symmetric_interval_hull(Phi2A_abs * symmetric_interval_hull(A * U))

    Ω0 = ConvexHull(X0, ϕ * X0 ⊕ δ*U ⊕ Eψ0 ⊕ Einit)
    Ud = δ*U ⊕ Eψ0
    In = IdentityMultiple(one(eltype(A)), size(A, 1))
    S_discr = ConstrainedLinearControlDiscreteSystem(ϕ, In, X, Ud)

    return InitialValueProblem(S_discr, Ω0)
end

# =================
# Correction hull
# =================

using IntervalMatrices: correction_hull

# TODO: add to IntervalMatrices
_correction_hull(A::IntervalMatrix, t, p) = correction_hull(A, t, p)
_correction_hull(A::AbstractMatrix, t, p) = correction_hull(IntervalMatrix(A), t, p)

# homogeneous case x' = Ax, x in X
# implements: Ω0 = CH(X0, exp(A*δ) * X0) ⊕ F*X0
# where F is the correction (interval) matrix
# if A is an interval matix, the exponential is overapproximated
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ::Float64, alg::CorrectionHullApproximation)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)
    @unpack order, exp_method = alg

    X0z = _convert_or_overapproximate(Zonotope, X0)
    if A isa IntervalMatrix
        Φ = exp_overapproximation(A, δ, order)
        Y = overapproximate(Φ * X0z, Zonotope)
    else
        Φ = _exp(A, δ, exp_method)
        Y = linear_map(Φ, X0z)
    end

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = _correction_hull(A, δ, order)
    R = overapproximate(F*X0z, Zonotope)
    Ω0 = minkowski_sum(H, R)

    ivp_discr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(ivp_discr, Ω0)
end

#=
function discretization(P::IVP{<:LCS, <:LazySet}, δ)
    # transforms to a CLCS (ConstrainedLinearContinuousSystem),
    # where the constraint is the universal set
    Snorm = Reachability.normalize(P.s)
    Pnorm = InitialValueProblem(Snorm, P.x0)
    return _discretize_homog(Pnorm, δ)
end
=#

#==
function discretization(P::IVP{<:CLCCS, <:LazySet}, δ)
    # transforms to a Constrained Linear Control Continuous System
    # in particular, inputs passed as a LazySet are wrapped into a
    # ConstantInput and inputs passed as a vector are wrapped as a
    # VaryingInput
    Snorm = Reachability.normalize(P.s)
    Pnorm = InitialValueProblem(Snorm, P.x0)
    return _discretize_inhomog(Pnorm, δ)
end
=#

# ======================

#=
function discretize(ivp_norm::IVP{<:CLCCS, <:LazySet}, δ::Float64, alg::ForwardApproximation)
    error("to-do")
end
=#

#=
function discretize(S::AbstractContinuousSystem, X0::LazySet, δ::Float64,
                    algo::AbstractApproximationModel=_default_approximation_model(ivp))
    # ...
    error("TODO")
end
=#
