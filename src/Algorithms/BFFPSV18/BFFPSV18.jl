"""
    BFFPSV18{N, ST, AM, IDX, BLK, RBLK, CBLK, PT, SP} <: AbstractContinuousPost

Implementation of the reachability method for linear systems using block decompositions.

## Fields

- `δ`             -- step-size of the discretization
- `setrep`        -- the set representation used in the decomposition
- `approx_model`  -- (optional, default: `Forward`) approximation model;
                     see `Notes` below for possible options
- `vars`          -- variables of interest
- `block_indices` -- tuple of integers to index each block that contains a variable of interest
- `row_blocks`    -- tuple of integer tuples to index variables associated to blocks of interest
- `column_blocks` -- tuple of integer tuples to index variables in the partition
- `partition`     -- tuple of integers that defines the decomposition of the
                     variable set `1, .., n` into given block sizes
- `sparse`        -- (optional, default: `false`) if `true`, assume that the state transition
                     matrix is sparse

See the `Examples` section below for some concrete examples of these options.

## Notes

This algorithm solves the set-based recurrence equation ``X_{k+1} = ΦX_k ⊕ V_k``
by using block decompositions. The algorithm was introduced in [[BFFPSV18]](@ref).

Comments about some fields:

- `N`    -- number type of the step-size, e.g. `Float64`
- `ST`   -- set representation used; this is either a concrete LazySet subtype,
            eg. `Interval{Float64, IntervalArithmetic.Interval{Float64}}`, or a tuple
            of concrete LazySet subtypes that is commensurate with the partition

The default approximation model is:

```julia
Forward(sih=:concrete, exp=:base, setops=:lazy)
```

### References

This algorithm is essentially an extension of the method in [[BFFPSV18]](@ref).
Blocks can have different dimensions and the set represenation can be different
for each block.

For a general introduction we refer to the dissertation [[SCHI18]](@ref).

Regarding the approximation model, by default we use an adaptation of the method
presented in [[FRE11]](@ref).
"""
struct BFFPSV18{N, ST, AM, IDX, BLK, RBLK, CBLK, PT, SP} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    vars::IDX
    block_indices::BLK
    row_blocks::RBLK
    column_blocks::CBLK
    partition::PT
    sparse::SP
end

# constructor with default approximation model
function BFFPSV18(; δ::N,
                    setrep::TST,
                    approx_model::AM=Forward(sih=:concrete, exp=:base, setops=:lazy),
                    vars::IDX,
                    block_indices::BLK,
                    row_blocks::RBLK,
                    column_blocks::CBLK,
                    partition::PT,
                    sparse::Bool=false
                 ) where {N, TST, AM, IDX, BLK, RBLK, CBLK, PT}

    # TODO add checks, alternative convenience constructors, etc.
    sparseval = Val(sparse)
    return BFFPSV18{N, setrep, AM, IDX, BLK, RBLK, CBLK, PT, typeof(sparseval)}(δ, approx_model,
                    vars, block_indices, row_blocks, column_blocks, partition, sparseval)
end

step_size(alg::BFFPSV18) = alg.δ
numtype(::BFFPSV18{N}) where {N} = N
setrep(::BFFPSV18{N, ST}) where {N, ST} = ST
rsetrep(alg::BFFPSV18{N, ST}) where {N, ST} = SparseReachSet{N, CartesianProductArray{N, ST}, length(alg.vars)}

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
