abstract type AbstractApproximationModel end
const AAModel = AbstractApproximationModel

export ForwardApproximation,
       BackwardApproximation,
       DiscreteApproximation

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
function _discretize_inhomog(ivp::IVP{<:CLCCS, <:LazySet}, δ::Float64, alg::ForwardApproximation)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    U = next_set(inputset(P), 1)
    ϕ = _exp(A, δ, alg.exp_method)
    A_abs = 
    Phi2Aabs = Reachability.ReachSets.Φ₂(abs.(A), δ, exp_method="base")
    Einit = symmetric_interval_hull(Phi2Aabs * symmetric_interval_hull((A * A) * X0))
    # cf. _discretize_interpolation_inhomog
    Eψ0 = symmetric_interval_hull(Phi2Aabs * symmetric_interval_hull(A * U0))
    Ω0 = ConvexHull(X0, ϕ * X0 ⊕ δ*U0 ⊕ Eψ0 ⊕ Einit)
    Ud = ConstantInput(δ*U0 ⊕ Eψ0)
    In = IdentityMultiple(one(Float64) * LinearAlgebra.I, size(A, 1))
    return IVP(CLCDS(ϕ, In, stateset(P.s), Ud), Ω0)
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

#const Partition{}

# concrete decomposition using a uniform block partition
#using LazySets.Arrays: projection_matrix

#const Partition{N, VT} = AbstractVector{VT} where {VT<:AbstractVector{Int}}

function _decompose(X::LazySet{N},
                    partition::AbstractVector{<:AbstractVector{Int}},
                    set_type::Type{ST}) where {N, ST<:LazySet}
    n = dim(X)
    result = Vector{ST}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        πS = Projection(S, block)
        result[i] = overapproximate(πS, X)
    end
    return CartesianProductArray(result)
end

# fallback implementation for conversion (if applicable) or overapproximation
function _convert_or_overapproximate(T::Type{<:AbstractPolytope}, X::LazySet)
    if applicable(convert, T, X)
        return convert(T, X)
    elseif applicable(overapproximate, X, T)
        return overapproximate(X, T)
    else
        return convert(T, overapproximate(X, Hyperrectangle))
    end
end

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
