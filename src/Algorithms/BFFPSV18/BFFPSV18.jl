"""
    BFFPSV18{N, ST, NI, IDX, BLK, RBLK CBLK} <: AbstractContinuousPost

Implementation of Girard - Le Guernic - Maler algorithm for reachability of
uncertain linear systems using zonotopes.

## Fields

- `δ`           -- step-size of the discretization
- `appro_model` -- (optional, default `_DEFAULT_APPROX_MODEL_GLGM06`) approximation model
                   for the discretization of the ODE; see `Notes` below
- `max_order`   -- (optional, default: `10`) maximum zonotope order

## Notes


## References

TODO: move these references to the general references

[1]
"""
@with_kw struct BFFPSV18{N, ST, AM, IDX, BLK, RBLK, CBLK, PT} <: AbstractContinuousPost
    δ::N
    setrep::ST  # remove ?
    approx_model::AM=_DEFAULT_APPROX_MODEL_BFFPSV18
    vars::IDX
    block_indices::BLK
    row_blocks::RBLK
    column_blocks::CBLK
    partition::PT
end

const _DEFAULT_APPROX_MODEL_BFFPSV18 = ForwardApproximation(sih_method="concrete", exp_method="base",
                                                            phi2_method="base") # TODO use Val{...}

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")

#=
opts1 = Options(:T => T,
               :δ => δ,
               :N => N,
               :vars=>[1,2,3,4,5],
               :partition => [[1], [2], [3], [4], [5]],
               #:set_type => Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}},
               :set_type => Interval{Float64, IA.Interval{Float64}},
               :row_blocks => [[1], [2], [3], [4], [5]],
               :column_blocks => [[1], [2], [3], [4], [5]],
               :num_type => Float64,
               :block_indices => [1, 2, 3, 4, 5],
               :sparse=>false);
=#
