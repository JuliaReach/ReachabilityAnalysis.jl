"""
    Disjunction <: Property

Type that represents a disjunction of properties.

### Fields

- `disjuncts` -- vector of properties (elements are reordered by this type)
- `reorder`   -- flag to indicate whether shuffling is allowed

### Notes

The following formula characterizes whether a set ``X`` satisfies a disjunction
``𝑃 = 𝑃_1 ∨ 𝑃_2 ∨ … ∨ 𝑃_m``:

```math
    X \\models 𝑃 \\iff X \\models 𝑃_j \\text{ for some } 1 ≤ j ≤ m
```

If the `reorder` flag is set, the disjuncts may be reordered after each call to
[`check`](@ref check(𝑃::Disjunction, X::LazySet{N}) where {N<:Real}) as a
heuristics to make subsequent checks faster.
"""
struct Disjunction <: Property
    disjuncts::Vector{Property}
    reorder::Bool
end

# default constructor with activated reordering
Disjunction(disjuncts::Vector{<:Property}) = Disjunction(disjuncts, true)

"""
    dim(𝑃::Disjunction)::Int

Return the dimension of a disjunction of properties.

### Input

- `𝑃` -- disjunction of properties

### Output

The dimension of the disjunction of properties.
"""
function dim(𝑃::Disjunction)::Int
    if isempty(𝑃.disjuncts)
        throw(ArgumentError("the dimension of an empty disjunction of " *
                            "properties is undefined"))
    end
    @inbounds return dim(𝑃.disjuncts[1])
end

"""
    check(𝑃::Disjunction, X::LazySet{N}; witness::Bool=false) where {N<:Real}

Check whether a convex set satisfies a disjunction of properties.

### Input

- `𝑃`       -- disjunction of properties
- `X`       -- convex set
- `witness` -- (optional, default: `false`) flag for returning a counterexample
               if the property is violated

### Output

* If `witness` option is deactivated: `true` iff `X` satisfies the property `𝑃`
* If `witness` option is activated:
  * `(true, [])` iff `X` satisfies the property `𝑃`
  * `(false, v)` iff `X` does not satisfy the property `𝑃` with witness `v`;
    note that `v == N[]` if 𝑃 is the empty disjunction

### Notes

By convention, the empty disjunction is equivalent to `false` and hence is
satisfied by no set.

If the `𝑃.reorder` flag is set, the disjuncts may be reordered as a heuristics
to make subsequent checks faster.
Since we check satisfaction from left to right, we move the disjunct for which
satisfaction was established to the front.

To be consistent with other propertes, the `witness` option only returns *one*
counterexample, namely for the left-most disjunct in the `disjuncts` vector.
We deactivate witness production for checking the remaining disjuncts.
"""
function check(𝑃::Disjunction, X::LazySet{N};
               witness::Bool=false) where {N<:Real}
    v = N[]
    create_witness = witness
    for (i, conjunct) in enumerate(𝑃.disjuncts)
        result = check(conjunct, X; witness=create_witness)
        if (create_witness && result[1]) || result
            _reorder!(𝑃, i)
            return witness ? (true, N[]) : true
        elseif create_witness
            v = result[2]
            # deactivate witness production for remaining checks
            create_witness = false
        end
    end
    return witness ? (false, v) : false
end

function _reorder!(𝑃::Disjunction, i::Int)
    if !𝑃.reorder || i == 1
        return nothing
    end
    first = 𝑃.disjuncts[i]
    while i > 1
        𝑃.disjuncts[i] = 𝑃.disjuncts[i-1]
        i -= 1
    end
    𝑃.disjuncts[1] = first
    return nothing
end

@inline function project(𝑃::Disjunction, vars::AbstractVector{Int})
    return Disjunction([project(disjunct, vars) for disjunct in 𝑃.disjuncts])
end
