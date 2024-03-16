using ..Overapproximate: _overapproximate

# =====================================
# Main solve for hybrid systems
# =====================================

const DTMAX_SIM_DEFAULT = 1e-1  # default time step for simulation with discrete events

# TODO API to specify the hybrid-specific options as a "DiscretePost" struct
# TODO pass internal options (eg. first_mode_representative, throw_error, etc.)

function solve(ivp::IVP{<:AbstractHybridSystem}, args...;
               max_jumps=1000, # maximum number of discrete transitions
               intersection_method::AbstractIntersectionMethod=HRepIntersection(), # method to take the concrete intersection in discrete transitions
               clustering_method::AbstractClusteringMethod=BoxClustering(), # method to perform clustering of the sets that cross a guard
               check_invariant_initial_states=false, # apply a disjointness check wrt each mode's invariant in the distribution of initial sets
               intersect_invariant_initial_states=false, # take the concrete intersection wrt each mode's invariant when distributing the initial states
               intersection_source_invariant_method=FallbackIntersection(), # method to take the concrete intersection with the source invariant
               first_mode_representative=true, # assume that the first mode is representative of the other modes when checking that the dimension in each mode is consistent
               intersect_source_invariant=true, # take the concrete intersection of the flowpipe with the source invariant
               disjointness_method=NoEnclosure(), # method to compute disjointness
               fixpoint_check=true, # if true, stop the integration when a fix point is detected
               WLtype=WaitingList, # type of waiting list used
               kwargs...)

    # distribute the initial condition across the different locations
    ivp_distributed = _distribute(ivp, WLtype; intersection_method=intersection_method,
                                  check_invariant=check_invariant_initial_states,
                                  intersect_invariant=intersect_invariant_initial_states)
    waiting_list = initial_state(ivp_distributed)
    H = system(ivp_distributed)

    # dimensional checks
    _check_dim(ivp_distributed; first_mode_representative=first_mode_representative)

    # get time span (or the emptyset if NSTEPS was specified)
    time_span = _get_tspan(args...; kwargs...)
    time_span0 = time_span
    # get the continuous post or find a default one
    cpost = _get_cpost(ivp_distributed, args...; kwargs...)
    if isnothing(cpost)
        cpost = _default_cpost(ivp_distributed, time_span; kwargs...)
    end

    # preallocate output flowpipe
    N = numtype(cpost)
    RT = rsetrep(cpost)
    out = Vector{Flowpipe{N,RT,Vector{RT}}}()
    sizehint!(out, max_jumps + 1)

    # list of (set, loc) tuples which have already been processed
    explored_list = typeof(waiting_list)()

    # preallocate output flowpipe strictly contained in each source invariant
    if intersect_source_invariant
        STwl = setrep(waiting_list)
        RTwl = ReachSet{N,STwl}
        out_in_inv = Vector{Flowpipe{N,RTwl,Vector{RTwl}}}()
        sizehint!(out_in_inv, max_jumps + 1)
    end

    # elapsed time accumulators
    t0 = tstart(time_span)
    @assert t0 == zero(t0) # we assume that the initial time is zero
    T = tend(time_span) # time horizon

    # counter for the number of transitions: using `count_jumps <= max_jumps` as
    # stopping criterion ensures that no more elements are added to the waiting
    # list after `max_jumps` discrete jumps
    count_jumps = 0

    while !isempty(waiting_list)
        (Δt0, elem) = pop!(waiting_list)
        tprev = tstart(Δt0)
        push!(explored_list, elem)

        # compute reachable states by continuous evolution
        q = location(elem)
        X0 = state(elem)
        S = mode(H, q)
        ivp_current = IVP(S, X0)
        time_span = TimeInterval(t0, T - tprev) # TODO generalization for t0 ≠ 0.. T-tprev+t0 ?
        F = post(cpost, ivp_current, time_span; Δt0=Δt0, kwargs...)

        # assign location q to this flowpipe
        F.ext[:loc_id] = q
        push!(out, F)

        I⁻ = stateset(H, q)

        #=
        # NOTE: here we may take the concrete intersection with the source invariant
        # We have to loop over each reach-set in F and only store the intersection
        # (possibly overapproximated) with the source invariant I⁻
        Two ways: (i) compute concrete intersection between F and I⁻
        (ii) compute disjointness between F and (I⁻)^C and only when it is non-zero
        compute the concrete intersection.
        =#
        if intersect_source_invariant
            # assign location q to this flowpipe
            F_in_inv = Flowpipe(undef, STwl, 0)

            for Ri in F
                # TODO refactor with reconstruct
                aux = _intersection(Ri, I⁻, intersection_source_invariant_method)
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
            r = target(H, t)

            # find reach-sets that may take the jump
            jump_rset_idx = Vector{Int}()
            for (i, X) in enumerate(F)
                _is_intersection_empty(X, G, disjointness_method) && continue
                push!(jump_rset_idx, i)
            end

            # continue if there is no guard enabled for this transition
            isempty(jump_rset_idx) && continue

            # apply clustering method to those sets which intersect the guard
            Xc = cluster(F, jump_rset_idx, clustering_method)

            for Xci in Xc
                # compute reachable states by discrete evolution
                X = apply(discrete_post, Xci, intersection_method)

                # do not add empty sets; checking `isempty` generalizes
                isa(X, EmptySet) && continue # isempty(X)

                count_jumps += 1
                hit_max_jumps = count_jumps > max_jumps
                if hit_max_jumps
                    @warn "maximum number of jumps reached; try increasing `max_jumps`"
                    continue
                end

                Xr = StateInLocation(X, r)
                if fixpoint_check
                    # check if this location has already been explored;
                    # if it is not the case, add it to the waiting list
                    add_to_waiting_list = !(Xr ⊆ explored_list)
                else
                    # We only push the set Xci if it intersects with time_span
                    add_to_waiting_list = !IA.isdisjoint(tspan(Xci), time_span0)
                end

                if add_to_waiting_list
                    push!(waiting_list, tspan(Xci), Xr)
                end
            end
        end # for
    end # while

    # wrap the flowpipe and algorithm in a solution structure
    if intersect_source_invariant
        HFout = HybridFlowpipe(out_in_inv)
    else
        HFout = HybridFlowpipe(out)
    end

    got_ensemble = get(kwargs, :ensemble, false)
    if got_ensemble
        @requires DifferentialEquations
        ensemble_sol = _solve_ensemble(ivp, args...; kwargs...)
        dict = Dict{Symbol,Any}(:ensemble => ensemble_sol)
        sol = ReachSolution(HFout, cpost, dict)
    else
        sol = ReachSolution(HFout, cpost)
    end
end

# =====================================
# Utility functions
# =====================================

# use the first_mode_representative flag to only check dimension of the first
# element in the waiting list
function _check_dim(ivp::InitialValueProblem{<:HybridSystem,<:AbstractWaitingList};
                    throw_error::Bool=true,
                    first_mode_representative=false)
    H = system(ivp)
    if first_mode_representative
        X0_1 = state(first(initial_state(ivp)))
        S1 = mode(H, 1)
        success = _check_dim(S1, X0_1; throw_error=throw_error)
    else
        success = true
        for elem in initial_state(ivp)
            X0_i = state(elem)
            Si = mode(H, location(elem))
            if !_check_dim(Si, X0_i; throw_error=throw_error)
                success = false
                break
            end
        end
    end
    return success
end

_distribute(ivp::IVP, WLtype::Type{<:WaitingList}; kwargs...) = _distribute(ivp; kwargs...)
function _distribute(ivp::IVP, WLtype::Type{<:MixedWaitingList}; kwargs...)
    return _distribute_mixed(ivp; kwargs...)
end

"""
    _distribute(ivp::InitialValueProblem{HS, ST};
                intersection_method::AbstractIntersectionMethod=nothing,
                check_invariant=false,
                intersect_invariant=false,
                ) where {HS<:HybridSystem, ST<:AdmissibleSet}

Distribute the set of initial states to each mode of a hybrid system.

### Input

- `system`              -- an initial value problem wrapping a mathematical system (hybrid)
                           and a set of initial states
- `intersection_method` -- method to perform the flowpipe-guard intersections
- `check_invariant`     -- (optional, default: `false`) if `false`, only add those modes
                           for which the intersection of the initial state with the invariant is non-empty
- `intersect_invariant` -- (optional, default: `false`) if `false`, take the concrete intersection with the invariant
 (and possibly overapproximate to keep `WaitingList` concretely typed)

### Output

A new initial value problem with the same hybrid system but where the set of initial
states is the list of tuples `(state, X0)`, for each state in the hybrid system.
"""
function _distribute(ivp::InitialValueProblem{HS,ST};
                     intersection_method=nothing,
                     check_invariant=false,
                     intersect_invariant=false) where {HS<:HybridSystem,ST<:AdmissibleSet}
    H = system(ivp)
    X0 = initial_state(ivp)
    N = eltype(X0)

    # NOTE using a the WaitingList, the set representation should be the same
    # for all sets => we have to convert (or overapproximate) the initial set X0
    # to the set representation that will be used for the sets after each discrete jump
    # we may be able to refactor this option, / ad an in-place _distribute
    # as this option is only used for storing (eventually overapproximating)
    # the new waiting list elements
    if !isnothing(intersection_method)
        STwl = setrep(intersection_method)
        X0 = _overapproximate(X0, STwl)
    else
        STwl = ST
    end

    waiting_list = WaitingList{TimeInterval,STwl,Int}()

    if !check_invariant
        for loc in states(H)
            push!(waiting_list, StateInLocation(X0, loc))
        end

    elseif check_invariant && !intersect_invariant
        for loc in states(H)
            Sloc = mode(H, loc)
            Y = stateset(Sloc)
            if !_is_intersection_empty(X0, Y)
                push!(waiting_list, StateInLocation(X0, loc))
            end
        end

    else # check_invariant && intersect_invariant
        for loc in states(H)
            Sloc = mode(H, loc)
            Y = stateset(Sloc)
            if !_is_intersection_empty(X0, Y)
                X0cut = intersection(X0, Y)
                X0cut_oa = overapproximate(X0cut, ST)
                push!(waiting_list, StateInLocation(X0cut_oa, loc))
            end
        end
    end
    return InitialValueProblem(H, waiting_list)
end

# in this implementation we use mixed waiting lists
function _distribute_mixed(ivp::InitialValueProblem{HS,ST};
                           intersection_method=nothing, # not used
                           check_invariant=false,
                           intersect_invariant=false) where {HS<:HybridSystem,ST<:AdmissibleSet}
    H = system(ivp)
    X0 = initial_state(ivp)
    N = eltype(X0)

    waiting_list = MixedWaitingList{TimeInterval,Vector{<:AdmissibleSet}}()

    if !check_invariant
        for loc in states(H)
            push!(waiting_list, StateInLocation(X0, loc))
        end

    elseif check_invariant && !intersect_invariant
        for loc in states(H)
            Sloc = mode(H, loc)
            Y = stateset(Sloc)
            if !_is_intersection_empty(X0, Y)
                push!(waiting_list, StateInLocation(X0, loc))
            end
        end

    else # check_invariant && intersect_invariant
        for loc in states(H)
            Sloc = mode(H, loc)
            Y = stateset(Sloc)
            if !_is_intersection_empty(X0, Y)
                X0cut = intersection(X0, Y)
                push!(waiting_list, StateInLocation(X0cut, loc))
            end
        end
    end
    return InitialValueProblem(H, waiting_list)
end

# the initial states are passed as a vector-of-tuples, each tuple being of the form
# (loc, X) where loc is an integer that corresponds to the mode and X is a set
function _distribute(ivp::InitialValueProblem{<:HybridSystem,Vector{Tuple{Int,ST}}};
                     intersection_method=nothing,
                     check_invariant=false,
                     intersect_invariant=false) where {ST<:AdmissibleSet}
    H = system(ivp)
    X0vec = initial_state(ivp) #  distributed initial states

    if !isnothing(intersection_method)
        STwl = setrep(intersection_method)
        X0vec = [(X0i[1], _overapproximate(X0i[2], STwl)) for X0i in X0vec]
    else
        STwl = ST
    end

    WL = WaitingList{TimeInterval,STwl,Int,StateInLocation{STwl,Int}}

    if !check_invariant && !intersect_invariant
        waiting_list = convert(WL, X0vec)
    else
        error("not implemented")
    end

    return InitialValueProblem(H, waiting_list)
end

# "duck-typing" the initial states passed as a vector-of-tuples, each tuple being
# of the form (X, loc) where loc is an integer that corresponds to the mode and X is a set
function _distribute(ivp::InitialValueProblem{HS,Vector{Tuple{ST,Int}}};
                     intersection_method=nothing,
                     check_invariant=false,
                     intersect_invariant=false) where {HS<:HybridSystem,ST<:AdmissibleSet}
    H = system(ivp)
    X0vec = initial_state(ivp) #  distributed initial states

    if !isnothing(intersection_method)
        STwl = setrep(intersection_method)
        X0vec = [(_overapproximate(X0i[1], STwl), X0i[2]) for X0i in X0vec]
    else
        STwl = ST
    end

    WL = WaitingList{TimeInterval,STwl,Int,StateInLocation{STwl,Int}}

    if !check_invariant && !intersect_invariant
        waiting_list = convert(WL, X0vec)
    else
        error("not implemented")
    end

    return InitialValueProblem(H, waiting_list)
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
                     ST, QT<:Tuple{Int}, ST}, VQT<:AbstractVector{QT}}
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
function constrained_dimensions(HS::HybridSystems.HybridSystem)
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
