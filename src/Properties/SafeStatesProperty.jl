"""
    SafeStatesProperty{N<:Real} <: Property

Type that represents a safety property characterized by a set of safe states.
The property is satisfied by a given set of states ``X`` if ``X`` is fully
contained in the set of safe states.

### Fields

- `safe`    -- convex set representing the safe states
- `witness` -- witness point (empty vector if not set)

### Notes

The following formula characterizes whether a set ``X`` satisfies a safety
property characterized by a set of safe states 𝑃:

```math
    X \\models 𝑃 \\iff X ⊆ 𝑃.\\texttt{safe}
```
"""
mutable struct SafeStatesProperty{N<:Real} <: Property
    safe::LazySet
    witness::Vector{N}

    SafeStatesProperty{N}(safe::LazySet) where {N<:Real} = new(safe, N[])
end

# type-less convenience constructor
SafeStatesProperty(safe::LazySet{N}) where {N<:Real} =
    SafeStatesProperty{N}(safe)

"""
    dim(𝑃::SafeStatesProperty)::Int

Return the dimension of a property with safe states.

### Input

- `𝑃` -- safety property with safe states

### Output

The dimension of the safe states.
"""
function dim(𝑃::SafeStatesProperty)::Int
    return dim(𝑃.safe)
end

"""
    check(𝑃::SafeStatesProperty, X::LazySet; witness::Bool=false)

Checks whether a convex set is contained in the set of safe states.

### Input

- `𝑃`       -- safety property with safe states
- `X`       -- convex set
- `witness` -- (optional, default: `false`) flag for returning a counterexample
               if the property is violated

### Output

Let ``Y`` be the safe states represented by 𝑃.
* If `witness` option is deactivated: `true` iff ``X ⊆ Y``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ Y``
  * `(false, v)` iff ``X ⊈ Y`` and ``v ∈ X \\setminus Y``
"""
@inline function check(𝑃::SafeStatesProperty, X::LazySet; witness::Bool=false)
    return ⊆(X, 𝑃.safe, witness)
end

@inline function project(𝑃::SafeStatesProperty, vars::AbstractVector{Int})
    proj = project(𝑃.safe, vars)
    return SafeStatesProperty(proj)
end
