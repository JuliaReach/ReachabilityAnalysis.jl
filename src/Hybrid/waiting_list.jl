# =====================================
# Structs to work with waiting lists
# =====================================

"""
    StateInLocation{ST, M<:Integer}

Struct that associates a set with a hybrid automaton's location index,
usually an integer.
"""
struct StateInLocation{ST,M}
    X::ST     # set representation
    loc_id::M # discrete location index
end

state(s::StateInLocation) = s.X
location(s::StateInLocation) = s.loc_id

# convert to a struct with broader type
function Base.convert(::Type{StateInLocation{ST,M}}, s::StateInLocation{WT,M}) where {ST,WT<:ST,M}
    return StateInLocation{ST,M}(s.X, s.loc_id)
end

# ===============================================================
# Abstract interface
# ===============================================================

"""
    AbstractWaitingList

Abstract supertype for all waiting list types.
"""
abstract type AbstractWaitingList end

# iterator interface
# array(w) should return the array of StateInLocation elements
@inline Base.getindex(w::AbstractWaitingList, i::Int) = getindex(array(w), i)
@inline Base.getindex(w::AbstractWaitingList, i::Number) = getindex(array(w), convert(Int, i))
@inline Base.getindex(w::AbstractWaitingList, I::AbstractVector) = getindex(array(w), I)
@inline Base.isempty(w::AbstractWaitingList) = isempty(array(w))
@inline Base.iterate(w::AbstractWaitingList) = iterate(array(w))
@inline Base.iterate(w::AbstractWaitingList, state) = iterate(array(w), state)
@inline Base.length(w::AbstractWaitingList) = length(array(w))
@inline Base.size(w::AbstractWaitingList) = (length(array(w)),)
@inline Base.first(w::AbstractWaitingList) = getindex(w, 1)
@inline Base.last(w::AbstractWaitingList) = getindex(w, lastindex(w))
@inline Base.firstindex(w::AbstractWaitingList) = 1
@inline Base.lastindex(w::AbstractWaitingList) = length(array(w))
@inline Base.eachindex(w::AbstractWaitingList) = eachindex(array(w))

# ===============================================================
# Waiting list with all sets of the same type
# ===============================================================

"""
    WaitingList{TN, ST, M, QT<:StateInLocation{ST, M}} <: AbstractWaitingList

Iterable container representing a list of pairs `(set, mode)` of a hybrid automaton.

### Fields

- `times`  -- vector with a time interval associated to each state
- `array`  -- vector of `StateInLocation`

### Notes

A `WaitingList` is a list of pairs ``(set_i, loc_i)` for `i in 1..k` where
`times` is a vector with a time interval associated to each state.

This waiting list allows for a unique set representation (`ST`) for all elements
of the list.
"""
struct WaitingList{TN,ST,M,QT<:StateInLocation{ST,M}} <: AbstractWaitingList
    times::Vector{TN}
    array::Vector{QT}

    function WaitingList(times::Vector{TN},
                         array::Vector{QT}) where {TN,ST,M,QT<:StateInLocation{ST,M}}
        @assert length(times) == length(array) || throw(ArgumentError("the lengths of the time " *
                                                                      "stamps and the array of sets should match, but they are $(length(times)) and $(length(array)) respectively"))

        return new{TN,ST,M,QT}(times, array)
    end
end

# constructor of empty waiting list
WaitingList{TN,ST,M}() where {TN,ST,M} = WaitingList(Vector{TN}(), Vector{StateInLocation{ST,M}}())
function WaitingList{TN,ST,M,QT}() where {TN,ST,M,QT<:StateInLocation{ST,M}}
    return WaitingList(Vector{TN}(), Vector{QT}())
end

# constructor without time-stamp
function WaitingList(array::Vector{QT}) where {ST,M,QT<:StateInLocation{ST,M}}
    return WaitingList{TimeInterval}(array)
end
function WaitingList{TN}(array::Vector{QT}) where {TN,ST,M,QT<:StateInLocation{ST,M}}
    times = fill(zeroI, length(array))
    return WaitingList(times, array)
end

# getter functions
@inline array(w::WaitingList) = w.array
@inline times(w::WaitingList) = w.times
setrep(w::WaitingList{TN,ST}) where {TN,ST} = ST
locrep(w::WaitingList{TN,ST,M}) where {TN,ST,M} = M
tspan(w::WaitingList, i) = w.times[i]

# we let elem::ET be possibly different than QT, for waiting lists which mixed set representations
@inline function Base.push!(w::WaitingList{TN,ST,M,QT}, elem::ET) where {TN,ST,M,QT,ET}
    push!(w.times, zeroI)
    push!(w.array, elem)
    return w
end
@inline function Base.push!(w::WaitingList{TN,ST,M,QT}, Δt::TN, elem::ET) where {TN,ST,M,QT,ET}
    push!(w.times, Δt)
    push!(w.array, elem)
    return w
end

@inline Base.pop!(w::WaitingList) = (pop!(w.times), pop!(w.array))

# containment check
function Base.:⊆(s::StateInLocation, w::WaitingList)
    contained = false
    q = location(s)
    X = state(s)
    for elem in w
        if (location(elem) == q) && (X ⊆ state(elem))
            contained = true
            break
        end
    end
    return contained
end

# conversion from vector-of-tuples to waiting list
function Base.convert(::Type{TW},
                      Q::Vector{Tuple{M,ST}}) where {TW<:WaitingList,M<:Integer,ST<:AdmissibleSet}
    waiting_list = TW()
    for Qi in Q # no intersection chcks
        q, Xi = Qi
        push!(waiting_list, StateInLocation(Xi, q))
    end
    return waiting_list
end

# "duck-typing" conversion from vector-of-tuples to waiting list
function Base.convert(::Type{TW},
                      Q::Vector{Tuple{ST,M}}) where {TW<:WaitingList,ST<:AdmissibleSet,M<:Integer}
    waiting_list = TW()
    for Qi in Q # no intersection chcks
        Xi, q = Qi
        push!(waiting_list, StateInLocation(Xi, q))
    end
    return waiting_list
end

# ===============================================================
# Waiting list allowing storage of sets of different types
# ===============================================================

# non-strictly typed waiting list, useful for cases in which the set reprsentation
# to be added in the list is not known a priori
struct MixedWaitingList{TN,QT<:StateInLocation} <: AbstractWaitingList
    times::Vector{TN}
    array::Vector{QT}

    function MixedWaitingList(times::Vector{TN}, array::Vector{QT}) where {TN,QT<:StateInLocation}
        @assert length(times) == length(array) || throw(ArgumentError("the lengths of the time " *
                                                                      "stamps and the array of sets should match, but they are $(length(times)) and $(length(array)) respectively"))

        return new{TN,QT}(times, array)
    end
end

# constructor of empty waiting list
MixedWaitingList{TN,QT}() where {TN,QT} = MixedWaitingList(Vector{TN}(), Vector{StateInLocation}())
