# =====================================
# Main solve for hybrid systems
# =====================================

# TODO API to specify the hybrid-specific options as a "DiscretePost" struct
# TODO pass internal options (eg. first_mode_representative, throw_error, etc.)

function solve(ivp::IVP{<:AbstractHybridSystem}, args...;
               max_jumps=100, # maximum number of discrete transitions
               intersection_method::AbstractIntersectionMethod=HRepIntersection(), # method to perform take the concrete intersection in discrete transitions
               clustering_method::AbstractClusteringMethod=BoxClustering(), # method to perform clustering of the sets that cross a guard
               check_invariant=false, # check intersection with invariant in the distribution of initial sets
               first_mode_representative=true, # assume that the first mode is representative of the other modes when checking that the dimension in each mode is consistent
               intersect_source_invariant=true, # take the concrete intersection of the flowpipe with the source invariant
               kwargs...)

    # distribute the initial condition across the different locations
    ivp_distributed = _distribute(ivp, intersection_method; check_invariant=check_invariant)
    waiting_list = initial_state(ivp_distributed)
    H = system(ivp_distributed)

    # dimensional checks
    _check_dim(ivp_distributed, first_mode_representative=first_mode_representative)

    # get time span (or the emptyset if NSTEPS was specified)
    time_span = _get_tspan(args...; kwargs...)

    # get the continuous post or find a default one
    cpost = _get_cpost(ivp_distributed, args...; kwargs...)
    if cpost == nothing
        cpost = _default_cpost(ivp_distributed, time_span; kwargs...)
    end

    # preallocate output flowpipe
    N = numtype(cpost)
    RT = rsetrep(cpost)
    out = Vector{Flowpipe{N, RT, Vector{RT}}}()
    sizehint!(out, max_jumps+1)

    # list of (set, loc) tuples which have already been processed
    STwl = setrep(waiting_list)
    MW = locrep(waiting_list)
    explored_list = WaitingList{N, STwl, MW}()

    # preallocate output flowpipe strictly contained in each source invariant
    if intersect_source_invariant
        RTwl = ReachSet{N, STwl}
        out_in_inv = Vector{Flowpipe{N, RTwl, Vector{RTwl}}}()
        sizehint!(out_in_inv, max_jumps+1)
    end

    # elapsed time accumulators
    t0 = tstart(time_span)
    @assert t0 == zero(t0) # NOTE: we assume that the initial time is zero
    T = tend(time_span) # time horizon

    count_jumps = 0

    while !isempty(waiting_list)  # && count_jumps <= max_jumps .. new stopping criterion?
        (tprev, elem) = pop!(waiting_list)
        push!(explored_list, elem)

        # compute reachable states by continuous evolution
        q = location(elem)
        X0 = state(elem)
        S = mode(H, q)
        time_span = TimeInterval(t0, T-tprev) # TODO generalization for t0 ≠ 0.. T-tprev+t0 ?
        F = post(cpost, IVP(S, X0), time_span; time_shift=tprev, kwargs...)

        # assign location q to this flowpipe
        F.ext[:loc_id] = q
        push!(out, F)

        #=
        # NOTE: here we may take the concrete intersection with the source invariant
        # We have to loop over each reach-set in F and only store the intersection
        # (possibly overapproximated) with the source invariant I⁻
        Two ways: (i) compute concrete intersection between F and I⁻
        (ii) compute disjointness between F and (I⁻)^C and only when it is non-zero
        compute the concrete intersection.
        =#
        if intersect_source_invariant
            I⁻ = stateset(H, q)
            # assign location q to this flowpipe
            F_in_inv = Flowpipe(undef, STwl, 0)
            for Ri in F
                # TODO refactor / reconstruct
                aux = intersection(set(Ri), I⁻)
                Raux = ReachSet(aux, tspan(Ri))
                push!(F_in_inv, Raux)
            end
            F_in_inv.ext[:loc_id] = q
            push!(out_in_inv, F_in_inv)
        end

        # process jumps for all outgoing transitions with source q
        for t in out_transitions(H, q)

            # instantiate post_d : X -> (R(X ∩ G ∩ I⁻) ⊕ W) ∩ I⁺
            discrete_post = DiscreteTransition(H, t)
            G = guard(discrete_post)

            # find reach-sets that may take the jump
            jump_rset_idx = Vector{Int}()
            for (i, X) in enumerate(F)
                _is_intersection_empty(X, G) && continue
                push!(jump_rset_idx, i)
            end

            # continue if there is no guard enabled for this transition
            isempty(jump_rset_idx) && continue

            # apply clustering method to those sets which intersect the guard
            # NOTE we assume that Clustering is applied and only Xc is only one reachset
            Xc = cluster(F, jump_rset_idx, clustering_method)
            tprev = tstart(Xc)

            # compute reachable states by discrete evolution
            X = apply(discrete_post, Xc, intersection_method) # .. may use intersection_method from @unwrap discrete post
            count_jumps += 1

            # check if this location has already been explored;
            # if it is not the case, add it to the waiting list
            r = target(H, t)
            Xr = StateInLocation(X, r)
            if (count_jumps <= max_jumps) && !(Xr ⊆ explored_list)
                push!(waiting_list, tprev, Xr)
                # NOTE save jumps_rset_idx, Xc? It might be interesting to store (for plots) the concrete
                # intersection with the source invariant F
            end
        end # for
    end # while

    # wrap the flowpipe and algorithm in a solution structure
    if intersect_source_invariant
        HFout = HybridFlowpipe(out_in_inv)
    else
        HFout = HybridFlowpipe(out)
    end
    sol = ReachSolution(HFout, cpost)
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
- `intersect_invariant` -- (optional, default: `false`) if `false`, take the concrete intersection with the invariant
 (and possibly overapproximate to keep `WaitingList` concretely typed)

### Output

A new initial value problem with the same hybrid system but where the set of initial
states is the list of tuples `(state, X0)`, for each state in the hybrid system.
"""
function _distribute(ivp::InitialValueProblem{HS, ST},
                     intersection_method::AbstractIntersectionMethod;
                     check_invariant=false,
                     intersect_invariant=false,
                     ) where {HS<:HybridSystem, N, ST<:LazySet{N}}
    H = system(ivp)
    X0 = initial_state(ivp)

    # NOTE using a the WaitingList, the set representation should be the same
    # for all sets => we have to convert (or overapproximate) the initial set X0
    # to the set representation that will be used for the sets after each discrete jump
    STwl = setrep(intersection_method)
    X0 = _overapproximate(X0, STwl)

    initial_states = WaitingList{N, STwl, Int}()

    if !check_invariant
        for loc in states(H)
            push!(initial_states, StateInLocation(X0, loc))
        end

    elseif check_invariant && !intersect_invariant
        for loc in states(H)
            Sloc = mode(H, loc)
            Y = stateset(Sloc)
            if !_is_intersection_empty(X0, Y)
                push!(initial_states, StateInLocation(X0, loc))
            end
        end

    else # check_invariant && intersect_invariant
        for loc in states(H)
            Sloc = mode(H, loc)
            Y = stateset(Sloc)
            if !_is_intersection_empty(X0, Y)
                X0cut = intersection(X0, Y)
                X0cut_oa = overapproximate(X0cut, ST)
                push!(initial_states, StateInLocation(X0cut_oa, loc))
            end
        end
    end
    return InitialValueProblem(H, initial_states)
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

# =====================================
# Discrete post structs
# =====================================

"""
    AbstractDiscretePost

Abstract supertype of all discrete post operators.
"""
abstract type AbstractDiscretePost <: AbstractPost end
