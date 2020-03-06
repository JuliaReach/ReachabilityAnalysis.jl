@with_kw struct BFFPSV18{#N<:Number,                 # numeric coefficient of the set rep
                         ST,            # parameter for the set type
                         NI, #<:Integer,               # integer type for the indices
                         IDX, #<:AbstractVector{NI},   # variable indices
                         BLK, #<:AbstractVector{NI},   # block indices
                         RBLK, #<:AbstractVector{NI},  # row indices
                         CBLK#<:AbstractVector{NI},  # column indices
                         } <: AbstractContinuousPost
    δ::Float64
    approx_model::AbstractApproximationModel=ForwardApproximation(sih_method="concrete",
                                                                  exp_method="base",
                                                                  phi2_method="base")

    #setrep::ST=Zonotope{Float64, Vector{Float64}, Matrix{Float64}}
    vars::IDX
    block_indices::BLK
    row_blocks::RBLK
    column_blocks::CBLK
    setrep::ST
    partition::AbstractVector{AbstractVector{Int}} # TODO add parameter
end

step_size(alg::BFFPSV18) = alg.δ
approx_model(alg::BFFPSV18) = alg.approx_model
setrep(alg::BFFPSV18{N, ST}) where {N, ST} = ST

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
