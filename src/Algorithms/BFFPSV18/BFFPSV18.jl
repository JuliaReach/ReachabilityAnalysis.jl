"""
    BFFPSV18{N, ST, AM, IDX, BLK, RBLK, CBLK, PT} <: AbstractContinuousPost

Implementation of reachability method for linear systems using block decompositions.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default: `Forward`) approximation model;
                    see `Notes` below for possible options
- `dim`          -- (optional default: `missing`) ambient dimension
- TODO: fix

## Notes

The type fields are:

- `N`         -- number type of the step-size
- `AM`        -- approximation model
- TODO: fix

The default approximation model is:

```julia
Forward(sih=:concrete, exp=:base, setops=:lazy)
```

This algorithm solves the set-based recurrence equation ``X_{k+1} = ΦX_k ⊕ V_k``
by using block decompositions. The algorithm was introduced in [[BFFPSV18]](@ref).

### References

This algorithm is essentially an extension of the method in [[BFFPSV18]](@ref).
Blocks can have different dimensions and the set represenation can be different
for each block.

For a general introduction we refer to the dissertation [[SCHI18]](@ref).

Regarding the approximation model, by default we use an adaptation of the method
presented in [[FRE11]](@ref).
"""
@with_kw struct BFFPSV18{N, ST, AM, IDX, BLK, RBLK, CBLK, PT} <: AbstractContinuousPost
    δ::N
    setrep::ST  # remove ?
    approx_model::AM=ForwardApproximation(sih=:concrete, exp=:base)
    vars::IDX
    block_indices::BLK
    row_blocks::RBLK
    column_blocks::CBLK
    partition::PT
end

step_size(alg::BFFPSV18) = alg.δ
numtype(::BFFPSV18{N}) where {N} = N

function rsetrep(::BFFPSV18{N}) where {N}
    error("not implemented")
end

#const _DEFAULT_APPROX_MODEL_BFFPSV18 =

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

#=
Example:

opts1 = Options(:T => T,
               :δ => δ,
               :N => N,
               :vars=>[1,2,3,4,5],
               :partition => [[1], [2], [3], [4], [5]],
               :setrep => Interval{Float64, IA.Interval{Float64}},
               :row_blocks => [[1], [2], [3], [4], [5]],
               :column_blocks => [[1], [2], [3], [4], [5]],
               :num_type => Float64,
               :block_indices => [1, 2, 3, 4, 5],
               :sparse=>false);
=#
