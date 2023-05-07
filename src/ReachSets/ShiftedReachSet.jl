# ================================================================
# Time-shifted reach-set
# ================================================================

struct ShiftedReachSet{N,RT<:AbstractLazyReachSet{N}} <: AbstractLazyReachSet{N}
    R::RT
    t0::N
end

# getter functions
@inline time_shift(srs::ShiftedReachSet) = srs.t0

# time domain interface
@inline tstart(srs::ShiftedReachSet) = tstart(srs.R) + time_shift(srs)
@inline tend(srs::ShiftedReachSet) = tend(srs.R) + time_shift(srs)
@inline tspan(srs::ShiftedReachSet) = TimeInterval(tstart(srs), tend(srs))
