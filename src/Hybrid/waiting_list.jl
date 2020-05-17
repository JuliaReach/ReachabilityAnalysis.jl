# =====================================
# Structs to work with waiting lists
# =====================================

# association of a set with a hybrid automaton location: (set, index of the mode)
struct StateInLocation{ST, M}
    X::ST # set representation
    loc_id::M # discrete location index
end

state(s::StateInLocation) = s.X
location(s::StateInLocation) = s.loc_id

# list of pairs (set_i, loc_i)
struct WaitingList{ST, M, QT<:StateInLocation{ST, M}}
    array::Vector{QT}
end

# constructor of empty waiting list
WaitingList{ST, M}() where {ST, M} = WaitingList(Vector{StateInLocation{ST, M}}())

# getter functions
@inline array(w::WaitingList) = w.array
setrep(w::WaitingList{ST}) where {ST} = ST
locrep(w::WaitingList{ST, M}) where {ST, M} = M

# iterator interface
@inline Base.getindex(w::WaitingList, i::Int) = getindex(array(w), i)
@inline Base.getindex(w::WaitingList, i::Number) = getindex(array(w), convert(Int, i))
@inline Base.getindex(w::WaitingList, I::AbstractVector) = getindex(array(w), I)
@inline Base.isempty(w::WaitingList) = isempty(array(w))
@inline Base.iterate(w::WaitingList) = iterate(array(w))
@inline Base.iterate(w::WaitingList, state) = iterate(array(w), state)
@inline Base.length(w::WaitingList) = length(array(w))
@inline Base.size(w::WaitingList) = (length(array(w)),)
@inline Base.first(w::WaitingList) = getindex(w, 1)
@inline Base.last(w::WaitingList) = getindex(w, lastindex(w))
@inline Base.firstindex(w::WaitingList) = 1
@inline Base.lastindex(w::WaitingList) = length(array(w))
@inline Base.eachindex(w::WaitingList) = eachindex(array(w))
@inline Base.push!(w::WaitingList, elem) = push!(w.array, elem)
@inline Base.pop!(w::WaitingList) = pop!(w.array)

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
