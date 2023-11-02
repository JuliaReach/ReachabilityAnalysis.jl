"""
    BFFPSV18{N, ST, AM, IDX, BLK, RBLK, CBLK} <: AbstractContinuousPost

Implementation of the reachability method for linear systems using block decompositions.

## Fields

- `δ`                -- step-size of the discretization
- `approx_model`     -- (optional, default: `Forward`) approximation model;
                        see `Notes` below for possible options
- `vars`             -- vector with the variables of interest
- `block_indices`    -- vector of integers to index each block that contains a variable of interest
- `row_blocks`       -- vector of integer vectors to index variables associated to blocks of interest
- `column_blocks`    -- vector of integer vectors to index variables in the partition
- `lazy_initial_set` -- (optional, default: `false`) if `true`, use a lazy decomposition of the initial states
                        after discretization
- `lazy_input`    -- (optional, default: `false`) if `true`, use a lazy decomposition of the input set
                      after discretization
- `sparse`        -- (optional, default: `false`) if `true`, assume that the state transition
                      matrix is sparse
- `view`          -- (optional, default: `false`) if `true`, use implementation that
                     uses arrays views

matrix is sparse

See the `Examples` section below for some concrete examples of these options.

## Notes

This algorithm solves the set-based recurrence equation ``X_{k+1} = ΦX_k ⊕ V_k``
by using block decompositions. The algorithm was introduced in [[BFFPSV18]](@ref).

Comments about some fields:

- `N`    -- number type of the step-size, e.g. `Float64`
- `ST`   -- set representation used; this is either a concrete LazySet subtype,
            eg. `Interval{Float64}`, or a tuple of concrete LazySet subtypes
            that is commensurate with the partition

The default approximation model is:

```julia
Forward(sih=:concrete, exp=:base, setops=:lazy)
```

TODO:

- clarify assumption about contiguous blocks

### Examples


### References

This algorithm is essentially an extension of the method in [[BFFPSV18]](@ref).
Blocks can have different dimensions and the set representation can be different
for each block.

For a general introduction we refer to the dissertation [[SCHI18]](@ref).

Regarding the approximation model, by default we use an adaptation of the method
presented in [[FRE11]](@ref).
"""
struct BFFPSV18{N,ST,AM,IDX,BLK,RBLK,CBLK} <: AbstractContinuousPost
    δ::N
    approx_model::AM
    vars::IDX
    block_indices::BLK
    row_blocks::RBLK
    column_blocks::CBLK
    lazy_initial_set::Bool
    lazy_input::Bool
    sparse::Bool
    view::Bool
end

function BFFPSV18(; δ::N,
                  setrep::ST=missing,
                  vars=missing,
                  partition::PT=missing,
                  dim::Union{Int,Missing}=missing,
                  approx_model::AM=Forward(; sih=:concrete, exp=:base, setops=:lazy),
                  lazy_initial_set::Bool=false,
                  lazy_input::Bool=false,
                  sparse::Bool=false,
                  view::Bool=true) where {N,ST,PT,AM}
    if !ismissing(dim)
        dimension = dim
    elseif !ismissing(partition)
        dimension = last(last(partition))
    else
        throw(ArgumentError("the constructor of `BFFPSV18` should specify at least the " *
                            "number of variables (`dim`) or the partition (`partition`)"))
    end

    if ismissing(vars) # if not specified, compute all variables
        vars = 1:dimension
    end

    if dimension == 1
        block_indices, row_blocks, column_blocks, setrep_tmp = _parse_opts_1D(vars, dimension,
                                                                              partition)
    else
        block_indices, row_blocks, column_blocks, setrep_tmp = _parse_opts_ND(vars, dimension,
                                                                              partition)
    end

    if ismissing(setrep)
        setrep = setrep_tmp
    end

    setrep = _concretize_setrep(setrep, N)

    return BFFPSV18{N,setrep,AM,typeof(vars),typeof(block_indices),typeof(row_blocks),
                    typeof(column_blocks)}(δ, approx_model,
                                           vars, block_indices, row_blocks, column_blocks,
                                           lazy_initial_set, lazy_input, sparse, view)
end

_concretize_setrep(setrep::Type{Interval}, N) = Interval{N}
_concretize_setrep(setrep::Type{Hyperrectangle}, N) = Hyperrectangle{N,Vector{N},Vector{N}}

# blocks of size 1
function _parse_opts_1D(vars, dim::Integer, partition)
    block_indices = collect(vars)
    row_blocks = [[i] for i in vars]
    column_blocks = [[i] for i in 1:dim]
    return block_indices, row_blocks, column_blocks, Interval
end

# blocks of size possibly greater than one
function _parse_opts_ND(vars, dim::Integer, partition)
    got_partition = !ismissing(partition)
    got_vars = !ismissing(vars)
    if got_partition
        column_blocks = copy(partition)
    elseif got_vars
        return _parse_opts_1D(vars, dim, partition)
    else
        column_blocks = [1:dim]
    end
    block_indices = Vector{Int}()
    row_blocks = Vector{Vector{Int}}()

    for v in vars
        for (j, bj) in enumerate(column_blocks)
            if v ∈ bj && j ∉ block_indices
                push!(block_indices, j)
                push!(row_blocks, bj)
                break
            end
        end
    end
    return block_indices, row_blocks, column_blocks, Hyperrectangle
end

# getter functions
step_size(alg::BFFPSV18) = alg.δ
numtype(::BFFPSV18{N}) where {N} = N
setrep(::BFFPSV18{N,ST}) where {N,ST} = ST

# the reduction is necessary to have the number of variables that are actually computed
function rsetrep(alg::BFFPSV18{N,ST}) where {N,ST}
    return SparseReachSet{N,CartesianProductArray{N,ST},length(reduce(vcat, alg.row_blocks))}
end

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
