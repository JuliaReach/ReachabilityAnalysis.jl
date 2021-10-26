# obs: S should normally be <:JuMP.MOI.AbstractOptimizer
struct ForwardBackward{EM, SO, SI, IT, BT, S} <: AbstractApproximationModel
    exp::EM
    setops::SO
    sih::SI
    inv::IT
    backend::BT
    solver::S
end

# backend for matrix exponential computations
hasbackend(alg::ForwardBackward) = !isnothing(alg.backend)
get_solver(alg::ForwardBackward) = alg.solver

# convenience constructor using symbols
function ForwardBackward(; exp=BaseExp, setops=:lazy, sih=:concrete, inv=false, backend=nothing, solver=nothing)
    isnothing(solver) && throw(ArgumentError("please specify an optimization solver, e.g. `solver=Ipopt.Optimizer()`"))
    return ForwardBackward(_alias(exp), _alias(setops), Val(sih), Val(inv), backend, solver)
end

function Base.show(io::IO, alg::ForwardBackward)
    print(io, "`ForwardBackward` approximation model with: \n")
    print(io, "    - exponentiation method: $(alg.exp) \n")
    print(io, "    - set operations method: $(alg.setops)\n")
    print(io, "    - symmetric interval hull method: $(alg.sih)\n")
    print(io, "    - invertibility assumption: $(alg.inv)")
    print(io, "    - polyhedral computations backend: $(alg.backend)")
end

Base.show(io::IO, m::MIME"text/plain", alg::ForwardBackward) = print(io, alg)

# ------------------------------------------------------------
# ForwardBackward Approximation: Homogeneous case
# ------------------------------------------------------------

function discretize(ivp::IVP{<:CLCS, <:LazySet}, δ, alg::ForwardBackward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    Φ = _exp(A, δ, alg.exp)

    A_abs = _elementwise_abs(A)
    Φcache = A == A_abs ? Φ : nothing
    P2A_abs = _Φ₂(A_abs, δ, alg.exp, alg.inv, Φcache)

    A² = A * A
    E₊ = sih(P2A_abs * sih(A² * X0, alg.sih), alg.sih)
    E₋ = sih(P2A_abs * sih((A² * Φ) * X0, alg.sih), alg.sih)

    n = size(A, 1)
    Eψ = ZeroSet(n)
    U = ZeroSet(n)
    Φᵀ = transpose(Φ)
    Ω0 = ContCH(δ, Φᵀ, X0, U, E₊, E₋, Eψ, get_solver(alg))

    # post-processing the lazy set
    Ω0 = _apply_setops(Ω0, alg)

    X = stateset(ivp)
    Sdis = ConstrainedLinearDiscreteSystem(Φ, X)
    return InitialValueProblem(Sdis, Ω0)
end

# ------------------------------------------------------------
# Forward Approximation: Inhomogeneous case
# ------------------------------------------------------------

function discretize(ivp::IVP{<:CLCCS, <:LazySet}, δ, alg::ForwardBackward)
    A = state_matrix(ivp)
    X0 = initial_state(ivp)
    Φ = _exp(A, δ, alg.exp)
    U = next_set(inputset(ivp), 1)

    A_abs = _elementwise_abs(A)
    Φcache = A == A_abs ? Φ : nothing
    P2A_abs = _Φ₂(A_abs, δ, alg.exp, alg.inv, Φcache)

    A² = A * A
    E₊ = sih(P2A_abs * sih(A² * X0, alg.sih), alg.sih)
    E₋ = sih(P2A_abs * sih((A² * Φ) * X0, alg.sih), alg.sih)
    Eψ = sih(P2A_abs * sih(A * U, alg.sih), alg.sih)

    Ud = δ*U ⊕ Eψ
    Φᵀ = transpose(Φ)
    Ω0 = ContCH(δ, Φᵀ, X0, U, E₊, E₋, Eψ, get_solver(alg))

    # post-processing the lazy set
    Ω0 = _apply_setops(Ω0, alg)

    X = stateset(ivp)
    In = IdentityMultiple(one(eltype(A)), size(A, 1))
    Sdis = ConstrainedLinearControlDiscreteSystem(Φ, In, X, Ud)
    return InitialValueProblem(Sdis, Ω0)
end

# struct to hold the set obtained with ForwardBackward discretization
# obs: OT is normally a MOI.AbstractOptimizer
mutable struct ContCH{N, MT, ST, UT, EPT, EMT, EST, OT} <: LazySet{N}
    δ::N
    Φᵀ::MT
    X0::ST
    U::UT
    E₊::EPT
    E₋::EMT
    Eψ::EST
    solver::OT
end

# constructor without solver specified
function ContCH(δ, Φᵀ, X0, U, E₊, E₋, Eψ; solver=nothing)
    return ContCH(δ, Φᵀ, X0, U, E₊, E₋, Eψ, solver)
end

function LazySets.dim(X::ContCH)
    return dim(X.X0)
end

# helper functions for the optimization solver
has_solver(X::ContCH) = !isnothing(X.solver)
get_solver(X::ContCH) = X.solver
function set_solver!(X::ContCH, solver)
    X.solver = solver
    return X
end

function _get_e(d, E)
    ee = σ(d, E)

    # revert sign whenever the direction is negative
    ee[d .< 0] .*= -1.
    return ee
end

# homogeneous case
function ω(λ, d, Φᵀ, X0, U::ZeroSet, Eψ, δ, e⁺, e⁻)
    aux1h = (1 - λ) * ρ(d, X0) + λ * ρ(Φᵀ*d, X0)
    aux2 = sum(min(λ * e⁺[i], (1 - λ)*e⁻[i]) * abs(d[i]) for i in eachindex(d))
    return aux1h + aux2
end

# inhomogeneous case
function ω(λ, d, Φᵀ, X0, U, Eψ, δ, e⁺, e⁻)
    aux1h = (1 - λ) * ρ(d, X0) + λ * ρ(Φᵀ*d, X0)
    aux1nh = λ * δ * ρ(d, U) + λ^2 * ρ(d, Eψ)
    aux2 = sum(min(λ * e⁺[i], (1 - λ)*e⁻[i]) * abs(d[i]) for i in eachindex(d))
    return aux1h + aux1nh + aux2
end

function load_forwardbackward_discretization()
return quote
    import .JuMP
    using .JuMP: MOI, Model, set_silent, register, optimize!, objective_value

function LazySets.ρ(d::AbstractVector, X::ContCH)
    @unpack δ, Φᵀ, X0, U, E₊, E₋, Eψ, solver = X

    !has_solver(X) && throw(ArgumentError("the optimization solver should be specified"))

    model = Model(typeof(get_solver(X)))
    set_silent(model)

    JuMP.@variable(model, 0 <= λ <= 1)

    e₊ = _get_e(d, E₊)
    e₋ = _get_e(d, E₋)

    _ω(λ) = ω(λ, d, Φᵀ, X0, U, Eψ, δ, e₊, e₋)
    register(model, :_ω, 1, _ω; autodiff = true)

    JuMP.@NLobjective(model, Max, _ω(λ))
    optimize!(model)
    return objective_value(model)
end

end # end quote
end # end load_forwardbackward_discretization
