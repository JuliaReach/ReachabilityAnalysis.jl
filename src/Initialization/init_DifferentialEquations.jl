using .DifferentialEquations
const DE = DifferentialEquations

using HybridSystems: states

# extend the solve API for initial-value problems

# =====================================
# Continuous system
# =====================================

function _solve_ensemble(ivp::InitialValueProblem, args...;
                         trajectories_alg=DE.Tsit5(),
                         ensemble_alg=DE.EnsembleThreads(),
                         inplace=true,
                         initial_states=nothing,
                         kwargs...)

    # get problem's vector field
    if inplace
        field = inplace_field!(ivp)
    else
        field = outofplace_field(ivp)
    end

    # get time span
    dt = _get_tspan(args...; kwargs...)
    tspan = (inf(dt), sup(dt))

    if isnothing(initial_states)
        # sample initial states
        X0 = initial_state(ivp)
        trajectories = get(kwargs, :trajectories, 100)
        initial_states = _sample_initial(X0, trajectories; kwargs...)
        # number of trajectories may increase if vertices got included
        trajectories = length(initial_states)
    else
        trajectories = length(initial_states)
    end

    # formulate ensemble ODE problem
    ensemble_prob = ODEProblem(field, first(initial_states), tspan)
    _prob_func(prob, i, repeat) = remake(prob, u0 = initial_states[i])

    # choose tolerances
    reltol = get(kwargs, :reltol, 1e-3)
    abstol = get(kwargs, :abstol, 1e-6)

    callback = get(kwargs, :callback, nothing)
    dtmax = get(kwargs, :dtmax, Inf)

    ensemble_prob = EnsembleProblem(ensemble_prob, prob_func = _prob_func)

    if isnothing(callback)
        result = DE.solve(ensemble_prob, trajectories_alg, ensemble_alg;
                          trajectories=trajectories, reltol=reltol,
                          abstol=abstol, dtmax=dtmax)
    else
        result = DE.solve(ensemble_prob, trajectories_alg, ensemble_alg;
                          trajectories=trajectories, reltol=reltol,
                          abstol=abstol, dtmax=dtmax, callback=callback)
    end
    return result
end

function _sample_initial(X0, trajectories; kwargs...)
    sampler = get(kwargs, :sampler, LazySets._default_sampler(X0))
    rng = get(kwargs, :rng, LazySets.GLOBAL_RNG)
    seed = get(kwargs, :seed, nothing)
    include_vertices = get(kwargs, :include_vertices, false)
    return sample(X0, trajectories; sampler=sampler, rng=rng,
                  seed=seed, include_vertices=include_vertices)
end

# =====================================
# Hybrid system
# =====================================

function _solve_ensemble(ivp::InitialValueProblem{<:AbstractHybridSystem},
                         args...; kwargs...)
    H = system(ivp)

    kwargs_sim = copy(kwargs)
    if :T in keys(kwargs_sim)
        # remove `T` because it is not allowed to define several time options
        pop!(kwargs_sim, :T)
    end

    time_span = _get_tspan(args...; kwargs...)
    t0_g = tstart(time_span)
    T = tend(time_span)

    termination_action = (integrator) -> terminate!(integrator)
    use_discrete_callback = get(kwargs, :use_discrete_callback, false)
    dtmax = get(kwargs, :dtmax, Inf)
    if use_discrete_callback && isinf(dtmax)
        # use a maximum time step when using discrete callbacks
        dtmax = DTMAX_SIM_DEFAULT
        kwargs_sim[:dtmax] = dtmax
    end

    jump_probability = get(kwargs, :jump_probability, 0.9)
    max_jumps = get(kwargs, :max_jumps, 100)

    # hybrid loop
    all_trajectories = []
    for initial_state in _sample_initial(ivp; kwargs...)
        loc, x0 = initial_state
        t0 = t0_g
        jump_count = 0
        trajectory_chain = []
        while true
            S_loc = mode(H, loc)
            ivp_loc = IVP(S_loc, x0)
            I⁻ = stateset(S_loc)

            # obtain trajectory in current mode
            if use_discrete_callback
                # use a discrete callback function based on membership
                condition = (u, t, integrator) -> u ∉ I⁻
                callback = DiscreteCallback(condition, termination_action)
            else
                # use a continuous callback function based on membership
                condition = (u, t, integrator) -> u ∉ I⁻ ? -1 : 1
                callback = ContinuousCallback(condition, termination_action)
            end
            trajectory = _solve_ensemble(ivp_loc, args...; initial_states=[x0],
                                         callback=callback, tspan=t0 .. T,
                                         kwargs_sim...)[1]

            # remove states from trajectory if they violate invariant
            # observation: this only concerns the last two steps
            if use_discrete_callback
                start_idx = max(1, length(trajectory.u) - 1)
                for idx in start_idx:length(trajectory.u)
                    if trajectory.u[idx] ∉ I⁻
                        indices_to_remove = idx:length(trajectory.u)
                        deleteat!(trajectory.u, indices_to_remove)
                        deleteat!(trajectory.t, indices_to_remove)
                        deleteat!(trajectory.k, indices_to_remove)
                        break
                    end
                end
            end

            # decide whether another jump should be attempted
            if jump_count > max_jumps
                # exceeded the allowed number of jumps
                try_jump = false
            elseif trajectory.t[end] >= T
                # time horizon reached → potentially do not take another jump
                try_jump = rand() <= jump_probability
            else
                try_jump = true
            end

            if try_jump
                # collect potential successors
                successors = []
                for t in out_transitions(H, loc)
                    G = guard(H, t)
                    asgn = resetmap(H, t)
                    loc′ = target(H, t)
                    I⁺ = stateset(H, loc′)

                    # find simulation steps that are enabled
                    # collect pairs (i, t, y) of indices i when transition t is
                    # enabled and leads to state y
                    for (i, xi) in enumerate(trajectory.u)
                        # check guard intersection
                        xi ∉ G && continue

                        # apply assignment
                        yi = MathematicalSystems.apply(asgn, xi)

                        # check target invariant
                        yi ∉ I⁺ && continue

                        # can take the transition
                        # store index, transition, and state
                        storage = (i, t, yi)
                        push!(successors, storage)
                    end
                end

                if isempty(successors)
                    jump_successful = false
                else
                    # choose a random successor
                    idx, t, y = rand(successors)
                    indices_to_remove = idx+1:length(trajectory.u)
                    deleteat!(trajectory.u, indices_to_remove)
                    deleteat!(trajectory.t, indices_to_remove)
                    deleteat!(trajectory.k, indices_to_remove)
                    loc′ = target(H, t)
                    push!(trajectory.u, y)
                    push!(trajectory.t, trajectory.t[end])
                    push!(trajectory.k, trajectory.k[end])
                    jump_successful = true
                end
            else
                jump_successful = false
            end

            # append new trajectory suffix to old trajectory prefix
            push!(trajectory_chain, trajectory)

            # test for termination
            if !jump_successful
                push!(all_trajectories, _merge_trajectory_chain(trajectory_chain))
                break
            end

            # continue previous simulation
            jump_count += 1
            loc = loc′
            x0 = trajectory.u[end]
            t0 = trajectory.t[end]
        end
    end
    return all_trajectories
end

# sample initial states of hybrid system given a list of pairs (loc, X0)
function _sample_initial(ivp::IVP{<:AbstractHybridSystem,
                                  <:Vector{<:Tuple{Integer, AdmissibleSet}}};
                         kwargs...)
    H = system(ivp)
    X0 = initial_state(ivp)

    # sample initial states from all possible initial regions
    trajectories = get(kwargs, :trajectories, 100)
    all_samples = []
    for (loc, X0_loc) in X0
        inv = stateset(mode(H, loc))
        X0_loc = X0_loc ∩ inv
        if isempty(X0_loc)
            continue
        end
        for x0 in _sample_initial(X0_loc, trajectories; kwargs...)
            push!(all_samples, (loc, x0))
        end
    end

    # filter for correct number of trajectories
    include_vertices = get(kwargs, :include_vertices, false)
    if !include_vertices
        rng = get(kwargs, :rng, LazySets.GLOBAL_RNG)
        samples = LazySets.Random.shuffle!(rng, all_samples)[1:trajectories]
    end

    return samples
end

# sample initial states of hybrid system given a set X0
function _sample_initial(ivp::IVP{<:AbstractHybridSystem, <:AdmissibleSet};
                         kwargs...)
    H = system(ivp)
    X0 = initial_state(ivp)
    X0_distributed = [(loc, X0) for loc in states(H)]
    return _sample_initial(IVP(H, X0_distributed))
end

# merge trajectory pieces into an ODESolution object
function _merge_trajectory_chain(trajectory_chain)
    prefix = trajectory_chain[1]
    @inbounds for i in 2:length(trajectory_chain)
        infix = trajectory_chain[i]
        append!(prefix.u, infix.u)
        append!(prefix.t, infix.t)
        append!(prefix.k, infix.k)
    end
    return prefix
end
