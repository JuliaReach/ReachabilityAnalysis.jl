"""
    inout_map_property(𝑃::PROPERTY,
                       partition::AbstractVector{<:AbstractVector{Int}},
                       blocks::AbstractVector{Int},
                       n::Int
                      )::PROPERTY where {PROPERTY<:Property}

Map a property to the dimensions of analyzed blocks.

### Input

- `𝑃`         -- property
- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition
- `n`         -- total number of input dimensions

### Output

A new property of reduced dimension, if necessary.

### Notes

The return type depends on the ambient dimension (`n`), the block indices in the
partition and on the type of property:

- The property `𝑃` is returned unchanged whenever `n` matches the dimensions
  in `blocks`.
- If the property is a `HalfSpace` (resp. `HPolyhedron`), the function `project`
  from `LazySets` returns a new property with a `HalfSpace` (resp. `HPolyhedron`)
  in the reduced dimensions, according to `blocks`.
- Otherwise, the dimensional reduction is implemented via a (lazy) `LinearMap`.
"""
function inout_map_property(𝑃::PROPERTY,
                            partition::AbstractVector{<:AbstractVector{Int}},
                            blocks::AbstractVector{Int},
                            n::Int
                           )::PROPERTY where {PROPERTY<:Property}

    # create a sorted list of all dimensions in `blocks` => available variables for projection
    proj = vcat(partition[blocks]...)

    # no change in the dimension => return the original property
    length(proj) == n && return 𝑃

    return inout_map_property_helper(𝑃, proj)
end

function inout_map_property_helper(𝑃::Conjunction, proj::Vector{Int})
    new_conjuncts = similar(𝑃.conjuncts)
    for (i, conjunct) in enumerate(𝑃.conjuncts)
        new_conjuncts[i] = inout_map_property_helper(conjunct, proj)
    end
    return Conjunction(new_conjuncts)
end

function inout_map_property_helper(𝑃::Disjunction, proj::Vector{Int})
    new_disjuncts = similar(𝑃.disjuncts)
    for (i, disjunct) in enumerate(𝑃.disjuncts)
        new_disjuncts[i] = inout_map_property_helper(disjunct, proj)
    end
    return Disjunction(new_disjuncts)
end

function inout_map_property_helper(𝑃::BadStatesProperty, proj::Vector{Int})
    return BadStatesProperty(project(𝑃.bad, proj))
end

function inout_map_property_helper(𝑃::SafeStatesProperty, proj::Vector{Int})
    return SafeStatesProperty(project(𝑃.safe, proj))
end

function compute_dimensions(partition, blocks)
    dimensions = Vector{Int}()
    dims = 0
    next_idx = 1
    for block_idx in blocks
        while next_idx < block_idx
            # sum up dimensions of skipped blocks
            dims += length(partition[next_idx])
            next_idx += 1
        end
        append!(dimensions, partition[next_idx])
        next_idx += 1
    end
    return dimensions
end
