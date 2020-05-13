# main solve for hybrid systems
function solve(ivp::IVP{<:AbstractHybridSystem}, args...;
               max_jumps=1000, kwargs...)


    # preliminary checks
    _check_dim(ivp)

    # get time span (or the emptyset if NSTEPS was specified)
    tspan = _get_tspan(args...; kwargs...)

    # get the continuous post or find a default one
    cpost = _get_cpost(ivp, args...; kwargs...)
    if cpost == nothing
        cpost = _default_cpost(ivp, tspan; kwargs...)
    end

    Q = _distribute(ivp)
    ST = setrep(cpost)

    waiting_list = WaitingList{N, ST}()

    while !isempty(waiting_list)
        # run the continuous-post operator
        F = post(cpost, ivp, tspan; kwargs...)
    end

    # wrap the flowpipe and algorithm in a solution structure
    sol = ReachSolution(F, cpost)
end

# ==========================================
# Convenience hybrid automaton constructors
# ==========================================

# hybrid automaton with one location and a self-loop
function HybridSystems.HybridSystem(mode::AbstractContinuousSystem, reset_map::AbstractMap)
    automaton = LightAutomaton(1)
    add_transition!(automaton, 1, 1, 1)
    return HybridSystem(automaton, [mode], [reset_map], [AutonomousSwitching()])
end

# ===================
# Utility functions
# ===================

# association of the pair (location, initial set)
struct StateInLocation{M, ST}
    loc_id::M # discrete location index
    X0::ST # set representation
end

# list of pairs (location, initial set)
struct WaitingList{M, ST, QT<:StateInLocation{M, ST}}
    array::Vector{QT}
end

# constructor of empty waiting list
WaitingList{M, ST}() where {M, ST} = WaitingList(Vector{StateInLocation{M,ST}}())

# getter functions
@inline array(wl::WaitingList) = wl.array

# iterator interface
@inline Base.getindex(wl::WaitingList, i::Int) = getindex(array(wl), i)
@inline Base.getindex(wl::WaitingList, i::Number) = getindex(array(wl), convert(Int, i))
@inline Base.getindex(wl::WaitingList, I::AbstractVector) = getindex(array(wl), I)
@inline Base.isempty(wl::WaitingList) = isempty(array(wl))
@inline Base.iterate(wl::WaitingList) = iterate(array(wl))
@inline Base.iterate(wl::WaitingList, state) = iterate(array(wl), state)
@inline Base.length(wl::WaitingList) = length(array(wl))
@inline Base.size(wl::WaitingList) = (length(array(wl)),)
@inline Base.first(wl::WaitingList) = getindex(wl, 1)
@inline Base.last(wl::WaitingList) = getindex(wl, lastindex(wl))
@inline Base.firstindex(wl::WaitingList) = 1
@inline Base.lastindex(wl::WaitingList) = length(array(wl))
@inline Base.eachindex(wl::WaitingList) = eachindex(array(wl))
@inline Base.push!(wl::WaitingList) = push!(wl.array)
@inline Base.pop!(wl::WaitingList) = pop!(wl.array)


"""
    _distribute(system::InitialValueProblem{<:HybridSystem, <:LazySet})

Distribute the set of initial states to each mode of a hybrid system.

### Input

- `system`          -- an initial value problem wrapping a mathematical system (hybrid)
                       and a set of initial states
- `check_invariant` -- (optional, default: `true`) if `true`, only add those modes
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
        initial_states = [(loc, X0) for loc in states(S)]
    else
        initial_states = Vector{Tuple{Int, ST}}()
        for loc in states(S)
            if !_is_intersection_empty(X0, stateset(loc))
                push!(initial_states, )
        end
    end
    return InitialValueProblem(S, initial_states)
end

# no-op if the initial set is given as tuples (location, set)
function _distribute(ivp::InitialValueProblem{HS, VQT};
                     check_invariant=false) where {HS<:HybridSystem,
                     ST, QT<:Tuple{<:Integer, ST}, VQT<:AbstractVector{QT}}
    !check_invariant && return ivp

end
