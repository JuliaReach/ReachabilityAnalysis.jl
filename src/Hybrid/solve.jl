# =====================================
# Main solve for hybrid systems
# =====================================

# TODO API to specify the hybrid-specific options as a "DiscretePost" struct
# TODO pass internal options (eg. first_mode_representative, throw_error, etc.)

function solve(ivp::IVP{<:AbstractHybridSystem}, args...;
               max_jumps=100,
               intersection::AbstractIntersectionMethod=HRepIntersection(),
               clustering::AbstractClusteringMethod=BoxClustering(),
               check_invariant=true,
               first_mode_representative=true,
               kwargs...)

    # distribute the initial condition across the different locations
    ivp_distributed = _distribute(ivp; check_invariant=check_invariant)
    waiting_list = initial_state(ivp_distributed)
    H = system(ivp_distributed)

    # dimensional checks
    _check_dim(ivp_distributed, first_mode_representative=first_mode_representative)

    # get time span (or the emptyset if NSTEPS was specified)
    tspan = _get_tspan(args...; kwargs...)

    # get the continuous post or find a default one
    cpost = _get_cpost(ivp_distributed, args...; kwargs...)
    if cpost == nothing
        cpost = _default_cpost(ivp_distributed, tspan; kwargs...)
    end

    # list of (set, loc) tuples which have already been processed
    STW = setrep(waiting_list)
    MW = locrep(waiting_list)
    explored_list = WaitingList{STW, MW}()

    # preallocate output flowpipe
    N = numtype(cpost)
    RT = rsetrep(cpost)
    out = Vector{Flowpipe{N, RT, Vector{RT}}}()
    sizehint!(out, max_jumps+1)

    # elapsed time accumulators
    t0 = tstart(tspan)
    #time_shift = zero(N)

    count_jumps = 0

    while !isempty(waiting_list) && count_jumps <= max_jumps
        elem = pop!(waiting_list)
        push!(explored_list, elem)

        # compute reachable states by continuous evolution
        q = location(elem)
        S = mode(H, q)
        X0 = state(elem)
        time_shift = 0.0
        F = post(cpost, IVP(S, X0), tspan; time_shift=time_shift, kwargs...)

        #time_shift += tend(F)  # update global time
        F.ext[:loc_id] = q      # assign location q to this flowpipe
        push!(out, F)

        #effective_tspan = tspan(F)

        # process jumps for all outgoing transitions with source q
        for t in out_transitions(H, q)

            # instantiate post_d : X -> (R(X ∩ G ∩ I⁻) ⊕ W) ∩ I⁺
            discrete_post = DiscreteTransition(H, t)
            G = guard(discrete_post)

            # find reach-sets that may take the jump
            jump_rset_idx = Vector{Int}()
            for (X, i) in enumerate(F)
                _is_intersection_empty(X, G) && continue
                push!(jump_rset_idx, i)
            end

            # continue if there is no guard enable for this transition
            isempty(jump_rset_idx) && continue

            # apply clustering method to those sets which intersect the guard
            Xc = cluster(F, jump_rset_idx, clustering)

            # compute reachable states by discrete evolution
            X = apply(discrete_post, Xc, method)

            # check if this location has already been explored;
            # if it is not the case, add it to the waiting list
            r = target(H, t)
            Xr = StateInLocation(X, r)
            if !(Xr ⊆ explored_list)
                push!(waiting_list, v)
            end
            count_jumps += 1
        end # for
    end # while

    # wrap the flowpipe and algorithm in a solution structure
    sol = ReachSolution(HybridFlowpipe(out), cpost)
end

# =====================================
# Utility functions
# =====================================

# use the first_mode_representative flag to only check dimension of the first
# element in the waiting list
function _check_dim(ivp::InitialValueProblem{<:HybridSystem, <:WaitingList};
                    throw_error::Bool=true,
                    first_mode_representative=false)

    H = system(ivp)
    if first_mode_representative
        X0_1 = state(first(initial_state(ivp)))
        S1 = mode(H, 1)
        success = _check_dim(S1, X0_1, throw_error=throw_error)
    else
        success = true
        for elem in initial_state(ivp)
            X0_i = state(elem)
            Si = mode(H, location(elem))
            if !_check_dim(Si, X0_i, throw_error=throw_error)
                success = false
                break
            end
        end
    end
    return success
end

"""
    _distribute(ivp::InitialValueProblem{HS, ST};
                check_invariant=false) where {HS<:HybridSystem, ST<:LazySet}

Distribute the set of initial states to each mode of a hybrid system.

### Input

- `system`          -- an initial value problem wrapping a mathematical system (hybrid)
                       and a set of initial states
- `check_invariant` -- (optional, default: `false`) if `false`, only add those modes
                       for which the intersection of the initial state with the invariant is non-empty

### Output

A new initial value problem with the same hybrid system but where the set of initial
states is the list of tuples `(state, X0)`, for each state in the hybrid system.
"""
function _distribute(ivp::InitialValueProblem{HS, ST};
                     check_invariant=false) where {HS<:HybridSystem, ST<:LazySet}
    S = system(ivp)
    X0 = initial_state(ivp)
    if !check_invariant
        initial_states = WaitingList([StateInLocation(X0, loc) for loc in states(S)])
    else
        initial_states = WaitingList{ST, Int}()
        for loc in states(S)
            Sloc = mode(S, loc)
            if !_is_intersection_empty(X0, stateset(Sloc))
                push!(initial_states, StateInLocation(X0, loc))
            end
        end
    end
    return InitialValueProblem(S, initial_states)
end

#=
function _distribute(ivp::InitialValueProblem{HS, WL},
                     check_invariant=false) where {HS<:HybridSystem, ST<:WaitingList}
    S = system(ivp)
    X0 = initial_state(ivp)
    if !check_invariant
        initial_states = WaitingList([StateInLocation(loc, X0) for loc in states(S)])
    else
        initial_states = WaitingList{ST, Int}()
        for loc in states(S)
            Sloc = mode(S, loc)
            if !_is_intersection_empty(X0, stateset(Sloc))
                push!(initial_states, StateInLocation(X0, loc))
            end
        end
    end
    return InitialValueProblem(S, initial_states)
end
=#

#=
# TODO distribute initial set exceptions:
# IA.Interval, IntervalBox, UninSet, UnionSetArray

# no-op if the ivp has a WaitingList

# initial set is given as a tuple (location, set)
# initial set is a vector of tuples [(loc_1, set_1), ..., (loc_n, set_n)]
# optionally it is checked that each set has a non-empty intersction with loc_i's invariant
function _distribute(ivp::InitialValueProblem{HS, VQT};
                     check_invariant=false) where {HS<:HybridSystem,
                     ST, QT<:Tuple{<:Integer, ST}, VQT<:AbstractVector{QT}}
    initial_states = WaitingList{Int, ST}()
    for loc in states(S)
        if check_invariant
            if !_is_intersection_empty(X0, stateset(loc))
                    push!(initial_states, StateInLocation(loc, X0))
                end
            end
        end
    end
end
=#

"""
    constrained_dimensions(HS::HybridSystem)::Dict{Int,Vector{Int}}

For each location, compute all dimensions that are constrained in the invariant
or the guard of any outgoing transition.

### Input

- `HS`  -- hybrid system

### Output

A dictionary mapping (`::Dict{Int,Vector{Int}}`) the index of each location
``ℓ`` to the dimension indices that are constrained in ``ℓ``.
"""
function constrained_dimensions(HS::HybridSystem)
    result = Dict{Int,Vector{Int}}()
    sizehint!(result, nstates(HS))
    for mode in states(HS)
        vars = Vector{Int}()
        append!(vars, constrained_dimensions(stateset(HS, mode)))
        for transition in out_transitions(HS, mode)
            append!(vars, constrained_dimensions(stateset(HS, transition)))
        end
        result[mode] = unique(vars)
    end

    return result
end
