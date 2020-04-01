using IntervalMatrices: correction_hull

"""
    AbstractApproximationModel

Abstract supertype for all approximation models.
"""
abstract type AbstractApproximationModel end

@with_kw struct Forward <: AbstractApproximationModel
    exp::Symbol=:base
    setops::Symbol=:lazy
    sih::Symbol=:concrete
    phi2::Symbol=:base
end

@with_kw struct Backward <: AbstractApproximationModel
    exp::Symbol=:base
    set::Symbol=:lazy
    sih::Symbol=:concrete
    phi2::Symbol=:base
end

# no bloating
struct Discrete <: AbstractApproximationModel
#
end

@with_kw struct CorrectionHull <: AbstractApproximationModel
   order::Int=10
   exp::Symbol=:base
end

function _default_approximation_model(ivp::IVP{<:AbstractContinuousSystem})
    return Forward()
end

# ============================================================
# Forward Approximation: Homogeneous case
# ============================================================

function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ::Float64, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    A_abs = _elementwise_abs(A)
    P2A_abs = Φ₂(A_abs, δ, alg.phi2)

    Ω0 = _discretize(A, X0, Φ, A_abs, P2A_abs, alg, Val(alg.setops))
    X = stateset(ivp)
    Sdiscr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdiscr, Ω0)
end

# change set_operations / setrep to -> ConvexHull (?)
function _discretize(A, X0, Φ, A_abs, P2A_abs, alg::Forward, setops::Val{:lazy})

    # "forward" algorithm, uses E⁺
    if alg.sih == :concrete
        # here the arguments to each symmetric_interval_hull are lazy
        Einit = symmetric_interval_hull(P2A_abs * symmetric_interval_hull((A * A) * X0))
    elseif alg.sih == :lazy
        Einit = SymmetricIntervalHull(P2A_abs * SymmetricIntervalHull((A * A) * X0))
    end

    Ω0 = ConvexHull(X0, Φ * X0 ⊕ Einit)
    return Ω0
end

function _discretize(A, X0, Φ, A_abs, P2A_abs, alg::Forward, setops::Val{:Interval})

    # "forward" algorithm, uses E⁺
    if alg.sih == :concrete
        # here the arguments to each symmetric_interval_hull are lazy
        Einit = symmetric_interval_hull(P2A_abs * symmetric_interval_hull((A * A) * X0))
    elseif alg.sih == :lazy
        Einit = SymmetricIntervalHull(P2A_abs * SymmetricIntervalHull((A * A) * X0))
    end

    Ω0 = ConvexHull(X0, Φ * X0 ⊕ Einit)
    return Ω0
end

function discretize(ivp::IVP{<:CLCS, Interval{N, IA.Interval{N}}}, δ::Float64, alg::Forward) where {N}
    A = state_matrix(ivp)
    @assert size(A, 1) == 1
    #@assert alg.setops == :Interval
    X0 = initial_state(ivp)

    a = A[1, 1]
    aδ = a * δ
    Φ = exp(aδ)
    A_abs = abs(a)

    # use inverse method
    @assert !iszero(a) "the given matrix should be invertible"

    # a_sqr = a * a
    #P2A_abs = (1/a_sqr) * (Φ - one(N) - aδ)
    #Einit = (P2A_abs * a_sqr) * RA._symmetric_interval_hull(X0).dat

    #P2A_abs = (1/a_sqr) * (Φ - one(N) - aδ)
    Einit = (Φ - one(N) - aδ) * _symmetric_interval_hull(X0).dat

    Ω0 = Interval(hull(X0.dat, Φ * X0.dat + Einit))
    X = stateset(ivp)
    # the system constructor creates a matrix
    Sdiscr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdiscr, Ω0)
end

# ============================================================
# Forward Approximation: Inhomogeneous case
# ============================================================

function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ::Float64, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)
    U = next_set(inputset(ivp), 1)
    ϕ = _exp(A, δ, alg.exp)
    A_abs = _elementwise_abs(A)
    Phi2A_abs = Φ₂(A_abs, δ, alg.phi2)

    @assert alg.sih == :concrete
    # TODO : specialize, add option to compute the concrete linear map
    Einit = symmetric_interval_hull(Phi2A_abs * symmetric_interval_hull((A * A) * X0))

    Eψ0 = symmetric_interval_hull(Phi2A_abs * symmetric_interval_hull(A * U))

    Ω0 = ConvexHull(X0, ϕ * X0 ⊕ δ*U ⊕ Eψ0 ⊕ Einit)
    Ud = δ*U ⊕ Eψ0
    In = IdentityMultiple(one(eltype(A)), size(A, 1))
    S_discr = ConstrainedLinearControlDiscreteSystem(ϕ, In, X, Ud)

    return InitialValueProblem(S_discr, Ω0)
end

# ===================================================
# Correction hull
# ===================================================

# homogeneous case x' = Ax, x in X
# implements: Ω0 = CH(X0, exp(A*δ) * X0) ⊕ F*X0
# where F is the correction (interval) matrix
# if A is an interval matix, the exponential is overapproximated
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ::Float64, alg::CorrectionHull)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)

    X0z = _convert_or_overapproximate(Zonotope, X0)
    if A isa IntervalMatrix
        Φ = exp_overapproximation(A, δ, alg.order)
        Y = overapproximate(Φ * X0z, Zonotope)
    else
        Φ = _exp(A, δ, alg.exp)
        Y = linear_map(Φ, X0z)
    end

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = correction_hull(A, δ, alg.order)
    R = overapproximate(F*X0z, Zonotope)
    Ω0 = _minkowski_sum(H, R)

    ivp_discr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(ivp_discr, Ω0)
end
