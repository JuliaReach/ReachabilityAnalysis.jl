"""
    BadStatesProperty{N<:Real} <: Property

Type that represents a safety property characterized by a set of bad states.
The property is satisfied by a given set of states if the intersection with the
set of bad states is empty.

### Fields

- `bad`     -- convex set representing the bad states
- `witness` -- witness point (empty vector if not set)

### Notes

The following formula characterizes whether a set ``X`` satisfies a safety
property characterized by a set of bad states 𝑃:

```math
    X \\models 𝑃 \\iff X ∩ 𝑃.\\texttt{bad} = ∅
```
"""
mutable struct BadStatesProperty{N<:Real} <: Property
    bad::LazySet
    witness::Vector{N}

    BadStatesProperty{N}(bad::LazySet) where {N<:Real} = new(bad, N[])
end

# type-less convenience constructor
BadStatesProperty(bad::LazySet{N}) where {N<:Real} =
    BadStatesProperty{N}(bad)

"""
    dim(𝑃::BadStatesProperty)::Int

Return the dimension of a property with bad states.

### Input

- `𝑃` -- safety property with bad states

### Output

The dimension of the bad states.
"""
function dim(𝑃::BadStatesProperty)::Int
    return dim(𝑃.bad)
end

"""
    check(𝑃::BadStatesProperty, X::LazySet; witness::Bool=false)

Checks whether a convex set is disjoint from the set of bad states.

### Input

- `𝑃`       -- safety property with bad states
- `X`       -- convex set
- `witness` -- (optional, default: `false`) flag for returning a counterexample
               if the property is violated

### Output

Let ``Y`` be the bad states represented by 𝑃.
* If `witness` option is deactivated: `true` iff ``X ∩ Y = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ Y = ∅``
  * `(false, v)` iff ``X ∩ Y ≠ ∅`` and ``v ∈ X ∩ Y``
"""
@inline function check(𝑃::BadStatesProperty, X::LazySet; witness::Bool=false)
    return isdisjoint(X, 𝑃.bad, witness)
end

@inline function project(𝑃::BadStatesProperty, vars::AbstractVector{Int})
    proj = project(𝑃.bad, vars)
    return BadStatesProperty(proj)
end
