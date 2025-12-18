# =======================================
# Flowpipe composition with a time-shift
# =======================================

"""
    ShiftedFlowpipe{FT<:AbstractFlowpipe, NT<:Number} <: AbstractFlowpipe

Type that lazily represents a flowpipe that has been shifted in time.

### Fields

- `F`  -- original flowpipe
- `t0` -- time shift

### Notes

This type can wrap any concrete subtype of `AbstractFlowpipe`, and the extra
field `t0` is such that the time spans of each reach-set in `F` are shifted
by the amount `t0` (which should be a subtype of `Number`).
"""
struct ShiftedFlowpipe{FT<:AbstractFlowpipe,NT<:Number} <: AbstractFlowpipe
    F::FT
    t0::NT
end

function ShiftedFlowpipe(vec::AbstractVector{<:AbstractLazyReachSet}, t0::Number)
    return ShiftedFlowpipe(Flowpipe(vec), t0)
end

# getter functions
@inline array(fp::ShiftedFlowpipe) = array(fp.F)
@inline flowpipe(fp::ShiftedFlowpipe) = fp.F
@inline time_shift(fp::ShiftedFlowpipe) = fp.t0

# time domain interface
@inline tstart(fp::ShiftedFlowpipe) = tstart(first(fp)) + fp.t0
@inline tend(fp::ShiftedFlowpipe) = tend(last(fp)) + fp.t0
@inline tspan(fp::ShiftedFlowpipe) = TimeIntervalC(tstart(fp), tend(fp))

@inline tstart(fp::ShiftedFlowpipe, i::Int) = tstart(fp.F[i]) + fp.t0
@inline tend(fp::ShiftedFlowpipe, i::Int) = tend(fp.F[i]) + fp.t0
@inline tspan(fp::ShiftedFlowpipe, i::Int) = TimeIntervalC(tstart(fp, i), tend(fp, i))

# TODO use interface?
project(fp::ShiftedFlowpipe, vars::AbstractVector) = project(fp, Tuple(vars))
project(fp::ShiftedFlowpipe; vars) = project(fp, Tuple(vars))

function project(fp::ShiftedFlowpipe, vars::NTuple{D,Int}) where {D}
    Xk = array(fp)
    # TODO: use projection of the reachsets
    if 0 ∈ vars # projection includes "time"
        # we shift the vars indices by one as we take the Cartesian prod with the time spans
        aux = vars .+ 1
        t0 = time_shift(fp)
        return map(X -> project(convert(Interval, tspan(X) + t0) × set(X), aux), Xk)
    else
        return map(X -> project(set(X), vars), Xk)
    end
end

# this method is analogue to project(::AbstractLazyReachSet, vars; check_vars=true)
# TODO add check_vars ?
function project(fp::ShiftedFlowpipe, i::Int, vars::NTuple{D,Int}) where {D}
    t0 = time_shift(fp)
    R = fp[i]
    if 0 ∈ vars
        # if the projection involves "time", we shift the vars indices by one as
        # we will take the Cartesian product of the reach-set with the time interval
        aux = vars .+ 1

        Δt = convert(Interval, tspan(R) + t0)
        proj = project(Δt × set(R), aux)
    else
        proj = project(set(R), vars)
    end

    return SparseReachSet(proj, tspan(R) + t0, vars)
end

# TODO: improve using mutable ShiftedFlowpipe struct (so that we can modify t0->t0+t1)
function shift(fp::ShiftedFlowpipe, t1::Number)
    return ShiftedFlowpipe(shift(fp.F, t1), fp.t0)
end
