"""
    Conjunction <: Property

Type that represents a conjunction of properties.

### Fields

- `conjuncts` -- vector of properties

### Notes

The following formula characterizes whether a set ``X`` satisfies a disjunction
``𝑃 = 𝑃_1 ∧ 𝑃_2 ∧ … ∧ 𝑃_m``:

```math
    X \\models 𝑃 \\iff X \\models 𝑃_j \\text{ for all } 1 ≤ j ≤ m
```
"""
struct Conjunction <: Property
    conjuncts::Vector{Property}
end

"""
    dim(𝑃::Conjunction)::Int

Return the dimension of a conjunction of properties.

### Input

- `𝑃` -- conjunction of properties

### Output

The dimension of the conjunction of properties.
"""
function dim(𝑃::Conjunction)::Int
    if isempty(𝑃.conjuncts)
        throw(ArgumentError("the dimension of an empty conjunction of " *
                            "properties is undefined"))
    end
    @inbounds return dim(𝑃.conjuncts[1])
end

"""
    check(𝑃::Conjunction, X::LazySet{N}; witness::Bool=false) where {N<:Real}

Check whether a convex set satisfies a conjunction of properties.

### Input

- `𝑃`       -- conjunction of properties
- `X`       -- convex set
- `witness` -- (optional, default: `false`) flag for returning a counterexample
               if the property is violated

### Output

* If `witness` option is deactivated: `true` iff `X` satisfies the property `𝑃`
* If `witness` option is activated:
  * `(true, [])` iff `X` satisfies the property `𝑃`
  * `(false, v)` iff `X` does not satisfy the property `𝑃` with witness `v`

### Notes

By convention, the empty conjunction is equivalent to `true` and hence is
satisfied by any set.
"""
function check(𝑃::Conjunction, X::LazySet{N};
               witness::Bool=false) where {N<:Real}
    for conjunct in 𝑃.conjuncts
        result = check(conjunct, X; witness=witness)
        if (witness && !result[1]) || !result
            return result
        end
    end
    return witness ? (true, N[]) : true
end

@inline function project(𝑃::Conjunction, vars::AbstractVector{Int})
    return Conjunction([project(conjunct, vars) for conjunct in 𝑃.conjuncts])
end
