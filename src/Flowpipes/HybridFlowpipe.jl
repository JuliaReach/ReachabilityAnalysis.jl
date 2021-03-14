# ================================
# Hybrid flowpipe
# ================================

"""
    HybridFlowpipe{N, D, FT<:AbstractFlowpipe, VOA<:VectorOfArray{N, D, Vector{FT}}} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of the same type, such that they are
contiguous in time.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes

The evaluation functions (in time) for this type do not assume that the flowpipes are contiguous in time.
That is, the final time of the `i`-th flowpipe does not match the start time of the `i+1`-th flowpipe.
"""
struct HybridFlowpipe{N, RT<:AbstractReachSet{N}, FT<:AbstractFlowpipe} <: AbstractFlowpipe
    Fk::VectorOfArray{RT, 2, Vector{FT}}
    ext::Dict{Symbol, Any}
end

function HybridFlowpipe(Fk::Vector{FT}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    ext = Dict{Symbol, Any}()
    return HybridFlowpipe{N, RT, FT}(voa, ext)
end

function HybridFlowpipe(Fk::Vector{FT}, ext::Dict{Symbol, Any}) where {N, RT<:AbstractReachSet{N}, FT<:Flowpipe{N, RT}}
    voa = VectorOfArray{RT, 2, Vector{FT}}(Fk)
    return HybridFlowpipe{N, RT, FT}(voa, ext)
end

function HybridFlowpipe(Fk::Vector{SFT}) where {N, RT, FT<:Flowpipe{N, RT}, NT<:Number, SFT<:ShiftedFlowpipe{FT, NT}}
    voa = VectorOfArray{RT, 2, Vector{SFT}}(Fk)
    ext = Dict{Symbol, Any}()
    return HybridFlowpipe{N, RT, SFT}(voa, ext)
end

# interface functions
array(fp::HybridFlowpipe) = fp.Fk
flowpipe(fp::HybridFlowpipe) = fp
numtype(::HybridFlowpipe{N}) where {N} = N
setrep(::Type{HybridFlowpipe{N, RT, FT}}) where {N, RT, FT} = setrep(RT)
setrep(::HybridFlowpipe{N, RT, FT}) where {N, RT, FT} = setrep(RT)
rsetrep(::Type{HybridFlowpipe{N, RT, FT}}) where {N, RT, FT} = RT
numrsets(fp::HybridFlowpipe) = mapreduce(length, +, fp)

# indexing: fp[j, i] returning the j-th reach-set of the i-th flowpipe
Base.getindex(fp::HybridFlowpipe, I::Int...) = getindex(fp.Fk, I...)

function tspan(fp::HybridFlowpipe)
    ti = minimum(tstart, fp)
    tf = maximum(tend, fp)
    return TimeInterval(ti, tf)
end

function Base.similar(fp::HybridFlowpipe{N, RT, FT}) where {N, RT, FT}
    return HybridFlowpipe(Vector{FT}())
end

function overapproximate(fp::HybridFlowpipe, args...)
    return HybridFlowpipe([overapproximate(F, args...) for F in fp], fp.ext)
end

function project(fp::HybridFlowpipe, args...)
    return [project(F, args...) for F in fp]
end

# LazySets interface
function LazySets.ρ(d::AbstractVector, fp::HybridFlowpipe)
    return maximum(ρ(d, F) for F in array(fp))
end

#function LazySets.σ(d::AbstractVector, fp::HybridFlowpipe)
#    error("not implemented")
#end

Base.:⊆(F::HybridFlowpipe, X::LazySet) = all(fp ⊆ X for fp in F)
Base.:⊆(F::HybridFlowpipe, Y::AbstractLazyReachSet) = all(fp ⊆ set(Y) for fp in F)

# evaluation for scalars
function (fp::HybridFlowpipe{N, RT})(t::Number) where {N, RT<:AbstractReachSet{N}}
    if t ∉ tspan(fp)
        throw(ArgumentError("time $t does not belong to the time span, " *
                            "$(tspan(fp)), of the given flowpipe"))
    end
    vec = Vector{RT}()
    for (k, fk) in enumerate(array(fp)) # loop over flowpipes
        if t ∈ tspan(fk)
            for Ri in fk  # loop over reach-sets for this flowpipe
                if t ∈ tspan(Ri)
                    push!(vec, Ri)
                end
            end
        end
    end
    return vec
end

# evaluation for time intervals
function (fp::HybridFlowpipe{N, RT})(dt::TimeInterval) where {N, RT<:AbstractReachSet{N}}
    if !(dt ⊆ tspan(fp)) # TODO IntervalArithmetic#409
        throw(ArgumentError("time interval $dt does not belong to the time span, " *
                            "$(tspan(fp)), of the given flowpipe"))
    end
    vec = Vector{RT}()
    for (k, fk) in enumerate(array(fp)) # loop over flowpipes
        if !isdisjoint(dt, tspan(fk))
            for R in fk(dt)  # loop over reach-sets for this flowpipe
                push!(vec, R)
            end
        end
    end
    return vec
end
