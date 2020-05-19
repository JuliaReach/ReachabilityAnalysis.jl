# =====================================
# Structs to work with waiting lists
# =====================================

# association of a set with a hybrid automaton location: (set, index of the mode)
struct StateInLocation{ST, M<:Integer}
    X::ST     # set representation
    loc_id::M # discrete location index
end

state(s::StateInLocation) = s.X
location(s::StateInLocation) = s.loc_id

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

# A WaitingList is a list of pairs (set_i, loc_i) for i in 1..k
# times is a vector with a time-stamp for each state
# This waiting list allows for a unique set representation (ST)
struct WaitingList{N, ST, M, QT<:StateInLocation{ST, M}} <: AbstractWaitingList
    times::Vector{N}
    array::Vector{QT}
    # TODO add inner constructor that checks that the lengths of times and array are the same

    function WaitingList(times::Vector{N}, array::Vector{QT}) where {N, ST, M, QT<:StateInLocation{ST, M}}
        @assert length(times) == length(array) || throw(ArgumentError("the lengths of the time " *
        "stamps and the array of sets should match, but they are $(length(times)) and $(length(array)) respectively"))

        return new{N, ST, M, QT}(times, array)
    end
end

# constructor of empty waiting list
WaitingList{N, ST, M}() where {N, ST, M} = WaitingList(Vector{N}(), Vector{StateInLocation{ST, M}}())

# constructor without time-stamp
WaitingList(array::Vector{QT}) where {ST, M, QT<:StateInLocation{ST, M}} = WaitingList{Float64}(array)
function WaitingList{N}(array::Vector{QT}) where {N, ST, M, QT<:StateInLocation{ST, M}}
    times = fill(zero(N), length(array))
    WaitingList(times, array)
end

# getter functions
@inline array(w::WaitingList) = w.array
@inline times(w::WaitingList) = w.times
setrep(w::WaitingList{N, ST}) where {N, ST} = ST
locrep(w::WaitingList{N, ST, M}) where {N, ST, M} = M
tstamp(w::WaitingList, i) = w.times[i]

@inline function Base.push!(w::WaitingList{N, ST, M, QT}, elem::QT) where {N, ST, M, QT}
    push!(w.times, zero(N))
    push!(w.array, elem)
    return w
end
@inline function Base.push!(w::WaitingList{N, ST, M, QT}, tstamp::N, elem::QT) where {N, ST, M, QT}
    push!(w.times, tstamp)
    push!(w.array, elem)
    return w
end
@inline Base.pop!(w::WaitingList) = (pop!(w.times), pop!(w.array))

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

#
# TODO MixedWaitingList : allow mixed types, such as small unions of
# StateInLocation with different set representations
#
