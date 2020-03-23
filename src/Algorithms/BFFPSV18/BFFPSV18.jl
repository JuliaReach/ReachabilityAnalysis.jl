"""
    BFFPSV18{N, ST, NI, IDX, BLK, RBLK CBLK} <: AbstractContinuousPost

Implementation of .... TODO

## Fields

TODO

## Notes

## References

TODO

"""
@with_kw struct BFFPSV18{N, ST, AM, IDX, BLK, RBLK, CBLK, PT} <: AbstractContinuousPost
    δ::N
    setrep::ST  # remove ?
    approx_model::AM=ForwardApproximation(sih_method=:concrete, exp_method=:base, phi2_method=:base)
    vars::IDX
    block_indices::BLK
    row_blocks::RBLK
    column_blocks::CBLK
    partition::PT
end

step_size(alg::BFFPSV18) = alg.δ

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
