const DEFAULT_TRAJECTORIES = 10

struct FixedStepSequence{N}
    F::Vector{Vector{N}}
    δ::N
end

# TODO plot recipe for FixedStepSequence

struct RungeKutta
    order::Int
    δ::Float64

    function RungeKutta(order=4, δ=1e-3)
        if order != 4
            throw(ArgumentError("currently, only order 4 is supported"))
        end
        return new(order, δ)
    end
end

function _solve_ensemble(ivp::IVP, args...;
                         trajectories_alg::RungeKutta=RungeKutta(4, 1e-3),
                         initial_states=nothing,
                         kwargs...)

    field = VectorField(ivp)

    # get time span
    dt = _get_tspan(args...; kwargs...)
    T = tend(dt) - tstart(dt)
    δ = trajectories_alg.δ
    r = range(zero(δ), stop=T, step=δ)
    order = Val(trajectories_alg.order)

    if isnothing(initial_states)
        # sample initial states
        X0 = initial_state(ivp)
        trajectories = get(kwargs, :trajectories, DEFAULT_TRAJECTORIES)
        initial_states = _sample_initial(X0, trajectories; kwargs...)
        # number of trajectories may increase if vertices got included
        trajectories = length(initial_states)
    else
        trajectories = length(initial_states)
    end

    N = Float64
    result = Vector{FixedStepSequence{N}}(undef, trajectories)
    @inbounds for i in 1:trajectories
        sequence = rk(field, initial_states[i], order, r)
        result[i] = FixedStepSequence(sequence, δ)
    end

    return result
end

function _sample_initial(X0, trajectories; kwargs...)
    sampler = get(kwargs, :sampler, LazySets._default_sampler(X0))
    rng = get(kwargs, :rng, GLOBAL_RNG)
    seed = get(kwargs, :seed, nothing)
    include_vertices = get(kwargs, :include_vertices, false)
    return sample(X0, trajectories; sampler=sampler, rng=rng,
                  seed=seed, include_vertices=include_vertices)
end

# Runge-Kutta method of order 4 with arbitrary but fixed steps
function rk(f, x0, order::Val{4}, steps::AbstractVector)
    x = x0
    results = Vector{typeof(x0)}(undef, length(steps))
    t0 = 0.0
    for (i, ti) in enumerate(steps)
        if ti > 0
            δ = ti - t0
            δ_half = δ / 2
            k1 = f(x)
            k2 = f(x + k1 * δ_half)
            k3 = f(x + k2 * δ_half)
            k4 = f(x + k3 * δ)
            x += 1/6 * δ * (k1 + 2 * k2 + 2 * k3 + k4)
            t0 = ti
        end
        results[i] = copy(x)
    end
    return results
end

# Runge-Kutta method with a step bound
function rk(f, x0, δ, order, k::Integer)
    r = range(zero(δ), step=δ, length=(k + 1))
    return rk(f, x0, order, r)
end

# Runge-Kutta method with a time bound
function rk(f, x0, δ, order, T::Real)
    r = range(zero(δ), stop=T, step=δ)
    return rk(f, x0, order, r)
end
