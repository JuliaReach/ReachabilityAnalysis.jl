struct ORBIT{N, VT, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM
end

step_size(alg::ORBIT) = alg.δ
numtype(::ORBIT{N}) where {N} = N

# convenience constructor using symbols
function ORBIT(; δ::N, approx_model::AM=NoBloating()) where {N, AM}
    VT = Vector{N}
    return ORBIT{N, VT, AM}(δ, approx_model)
end

function setrep(::ORBIT{N, VT, AM}) where {N, VT, AM}
    Singleton{N, VT}
end

function rsetrep(alg::ORBIT{N}) where {N}
    ST = sertrep(alg)
    ReachSet{N, ST}
end

include("post.jl")
