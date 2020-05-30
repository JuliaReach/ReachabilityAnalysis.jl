using IntervalMatrices: correction_hull

"""
    AbstractApproximationModel

Abstract supertype for all approximation models.
"""
abstract type AbstractApproximationModel end

# "forward" algorithm, uses E⁺
struct Forward{EM, SO, SI} <: AbstractApproximationModel
    exp::EM
    setops::SO
    sih::SI
end

# convenience constructor using symbols
function Forward(; exp::Symbol=:base, setops::Symbol=:lazy, sih::Symbol=:concrete)
    Forward(Val(exp), Val(setops), Val(sih))
end

# TODO: improve display of Forward()

struct Backward{EM, SO, SI} <: AbstractApproximationModel
    exp::EM
    setops::SO
    sih::SI
end

# convenience constructor using symbols
function Backward(; exp::Symbol=:base, setops::Symbol=:lazy, sih::Symbol=:concrete)
    Backward(Val(exp), Val(setops), Val(sih))
end

# no bloating or "discrete time" approximation model ref Eqs. (14) in [[BFFPSV18]]
struct NoBloating{EM, SO} <: AbstractApproximationModel
    exp::EM
    setops::SO
end

# convenience constructor using symbols
function NoBloating(; exp::Symbol=:base, setops::Symbol=:lazy)
    NoBloating(Val(exp), Val(setops))
end

struct CorrectionHull{EM} <: AbstractApproximationModel
   order::Int
   exp::EM
end

# convenience constructor using symbols
function CorrectionHull(; order::Int=10, exp::Symbol=:base)
    CorrectionHull(order, Val(exp))
end

function _default_approximation_model(ivp::IVP{<:AbstractContinuousSystem})
    return Forward()
end

# ============================================================
# Forward Approximation: Homogeneous case
# ============================================================

sih(X, ::Val{:lazy}) where {EM, SO} = SymmetricIntervalHull(X)
sih(X, ::Val{:concrete}) where {EM, SO} = _symmetric_interval_hull(X)

function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    A_abs = _elementwise_abs(A)
    P2A_abs = Φ₂(A_abs, δ, alg.exp)

    Einit = sih(P2A_abs * sih((A * A) * X0, alg.sih), alg.sih)
    Ω0 = ConvexHull(X0, Φ * X0 ⊕ Einit)
    X = stateset(ivp)
    Sdiscr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdiscr, Ω0)
end

#=
TODO remove
function _get_Einit(A, X0, P2A_abs, ::Forward{EM, Val{:lazy}, Val{:lazy}}) where {EM}
    Einit = SymmetricIntervalHull(P2A_abs * SymmetricIntervalHull((A * A) * X0))
end

function _get_Einit(A, X0, P2A_abs, ::Forward{EM, Val{:lazy}, Val{:concrete}}) where {EM}
    Einit = symmetric_interval_hull(P2A_abs * symmetric_interval_hull((A * A) * X0))
end

function _discretize(A, X0, P2A_abs, alg::Forward{EM, Val{:interval}, Val{:concrete}}) where {EM}
    Einit = symmetric_interval_hull(P2A_abs * symmetric_interval_hull((A * A) * X0))
end

function _discretize(A, X0, P2A_abs, alg::Forward{EM, Val{:interval}, Val{:lazy}}) where {EM}
    Einit = SymmetricIntervalHull(P2A_abs * SymmetricIntervalHull((A * A) * X0))
end
=#

function discretize(ivp::IVP{<:CLCS, Interval{N, IA.Interval{N}}}, δ, alg::Forward) where {N}
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

# TODO : specialize, add option to compute the concrete linear map
function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ, alg::Forward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    A_abs = _elementwise_abs(A)
    Phi2A_abs = Φ₂(A_abs, δ, alg.exp)

    Einit = sih(Phi2A_abs * sih((A * A) * X0, alg.sih), alg.sih)

    U = next_set(inputset(ivp), 1) # inputset(ivp)
    Eψ0 = sih(Phi2A_abs * sih(A * U, alg.sih), alg.sih)

    Ud = δ*U ⊕ Eψ0
    In = IdentityMultiple(one(eltype(A)), size(A, 1))

    Ω0 = ConvexHull(X0, Φ * X0 ⊕ Ud ⊕ Einit)
    X = stateset(ivp)
    Sdiscr = ConstrainedLinearControlDiscreteSystem(Φ, In, X, Ud)
    return InitialValueProblem(Sdiscr, Ω0)
end

# ===================================================
# Correction hull
# ===================================================

# homogeneous case x' = Ax, x in X
# implements: Ω0 = CH(X0, exp(A*δ) * X0) ⊕ F*X0
# where F is the correction (interval) matrix
# if A is an interval matix, the exponential is overapproximated
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::CorrectionHull)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    X = stateset(ivp)

    X0z = _convert_or_overapproximate(Zonotope, X0)
    if A isa IntervalMatrix
        Φ = exp_overapproximation(A, δ, alg.order)

        #Φ = IntervalMatrices.scale_and_square(A, 10, δ, 10)
        Y = _overapproximate(Φ * X0z, Zonotope)
    else
        Φ = _exp(A, δ, alg.exp)
        Y = _linear_map(Φ, X0z)
    end

    H = overapproximate(CH(X0z, Y), Zonotope)
    F = correction_hull(A, δ, alg.order)
    R = _overapproximate(F * X0z, Zonotope)
    Ω0 = _minkowski_sum(H, R)

    ivp_discr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(ivp_discr, Ω0)
end

# ============================================================
# NoBloating Approximation
# ============================================================

# homogeneous case
function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::NoBloating)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)

    Ω0 = copy(X0)
    X = stateset(ivp)
    Sdiscr = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdiscr, Ω0)
end

# inhomogeneous case
function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ, alg::NoBloating)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)

    Φ = _exp(A, δ, alg.exp)
    M = Φ₁(A, δ, alg.exp)
    U = next_set(inputset(ivp), 1) # inputset(ivp)
    Ud = M * U # TODO assert alg.setops  / alg.inputsetops are lazy
    In = IdentityMultiple(one(eltype(A)), size(A, 1))

    Ω0 = copy(X0)
    X = stateset(ivp)
    Sdiscr = ConstrainedLinearControlDiscreteSystem(Φ, In, X, Ud)
    return InitialValueProblem(Sdiscr, Ω0)
end
