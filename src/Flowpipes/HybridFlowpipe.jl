# ================================
# Hybrid flowpipe
# ================================

"""
    HybridFlowpipe{N,RT<:AbstractReachSet{N},FT<:AbstractFlowpipe} <: AbstractFlowpipe

Type that wraps a vector of flowpipes of the same type, such that they are
contiguous in time.

### Fields

- `Fk`  -- vector of flowpipes
- `ext` -- (optional, default: empty) dictionary for extensions

### Notes

The evaluation functions (in time) for this type do not assume that the flowpipes are contiguous in time.
That is, the final time of the `i`-th flowpipe does not match the start time of the `i+1`-th flowpipe.
"""
struct HybridFlowpipe{N,RT<:AbstractReachSet{N},FT<:AbstractFlowpipe} <: AbstractFlowpipe
    Fk::Vector{FT}
    ext::Dict{Symbol,Any}
end

function HybridFlowpipe(Fk::Vector{FT},
                        ext::Dict{Symbol,Any}=Dict{Symbol,Any}()) where {N,RT<:AbstractReachSet{N},
                                                                         FT<:Flowpipe{N,RT}}
    return HybridFlowpipe{N,RT,FT}(Fk, ext)
end

function HybridFlowpipe(Fk::Vector{SFT},
                        ext::Dict{Symbol,Any}=Dict{Symbol,Any}()) where {N,RT,FT<:Flowpipe{N,RT},
                                                                         NT<:Number,
                                                                         SFT<:ShiftedFlowpipe{FT,
                                                                                              NT}}
    return HybridFlowpipe{N,RT,SFT}(Fk, ext)
end

# interface functions
array(fp::HybridFlowpipe) = fp.Fk
flowpipe(fp::HybridFlowpipe) = fp
numtype(::HybridFlowpipe{N}) where {N} = N
setrep(::Type{HybridFlowpipe{N,RT,FT}}) where {N,RT,FT} = setrep(RT)
setrep(::HybridFlowpipe{N,RT,FT}) where {N,RT,FT} = setrep(RT)
rsetrep(::Type{HybridFlowpipe{N,RT,FT}}) where {N,RT,FT} = RT
numrsets(fp::HybridFlowpipe) = mapreduce(length, +, fp)

# indexing: fp[j, i] returning the j-th reach-set of the i-th flowpipe
Base.getindex(fp::HybridFlowpipe, I::Int...) = getindex(fp.Fk, I...)

function tspan(fp::HybridFlowpipe)
    ti = minimum(tstart, fp)
    tf = maximum(tend, fp)
    return TimeIntervalC(ti, tf)
end

function Base.similar(fp::HybridFlowpipe{N,RT,FT}) where {N,RT,FT}
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
function (fp::HybridFlowpipe{N,RT})(t::Number) where {N,RT<:AbstractReachSet{N}}
    if !IA.in_interval(t, tspan(fp))
        throw(ArgumentError("time $t does not belong to the time span, " *
                            "$(tspan(fp)), of the given flowpipe"))
    end
    vec = Vector{RT}()
    for (k, fk) in enumerate(array(fp)) # loop over flowpipes
        if IA.in_interval(t, tspan(fk))
            for Ri in fk  # loop over reach-sets for this flowpipe
                if IA.in_interval(t, tspan(Ri))
                    push!(vec, Ri)
                end
            end
        end
    end
    return vec
end

# evaluation for time intervals
function (fp::HybridFlowpipe{N,RT})(dt::TimeInterval) where {N,RT<:AbstractReachSet{N}}
    if !(dt ⊆ tspan(fp)) # ⊊ doesn't work see, IntervalArithmetic#409
        throw(ArgumentError("time interval $dt does not belong to the time span, " *
                            "$(tspan(fp)), of the given flowpipe"))
    end
    vec = Vector{RT}()
    for (k, fk) in enumerate(array(fp)) # loop over flowpipes
        dtint = IA.intersect_interval(dt, tspan(fk))
        if !IA.isempty_interval(dtint)
            for R in fk(dtint)  # loop over reach-sets for this flowpipe
                push!(vec, R)
            end
        end
    end
    return vec
end
