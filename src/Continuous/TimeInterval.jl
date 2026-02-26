struct TimeInterval{T<:Union{IA.Interval,Interval,EmptySet}}
    i::T
end

function TimeInterval(t0::Real, t1::Real)
    return TimeInterval(IA.interval(t0, t1))
end

const zeroT = TimeInterval(zeroI)

convert(::Type{Interval}, Δt::TimeInterval{<:Interval}) = Δt.i
convert(::Type{Interval}, Δt::TimeInterval{<:IA.Interval}) = Interval(Δt.i)

convert(::Type{IA.Interval}, Δt::TimeInterval{<:Interval}) = Δt.i.dat
convert(::Type{IA.Interval}, Δt::TimeInterval{<:IA.Interval}) = Δt.i

tstart(Δt::TimeInterval{<:Interval}) = low(Δt.i, 1)
tstart(Δt::TimeInterval{<:IA.Interval}) = inf(Δt.i)

tend(Δt::TimeInterval{<:Interval}) = high(Δt.i, 1)
tend(Δt::TimeInterval{<:IA.Interval}) = sup(Δt.i)

tspan(Δt::TimeInterval) = Δt

diam(Δt::TimeInterval{<:Interval}) = diameter(Δt.i)
diam(Δt::TimeInterval{<:IA.Interval}) = IA.diam(Δt.i)

isempty(Δt::TimeInterval) = isempty(Δt.i)

in(x::Real, Δt::TimeInterval) = in(x, Δt.i)

function _isapprox(Δt::TimeInterval, Δs::TimeInterval)
    return (tstart(Δt) ≈ tstart(Δs)) && (tend(Δt) ≈ tend(Δs))
end

@commutative function +(Δt::TimeInterval{<:Interval}, t::Real)
    return TimeInterval(Interval(tstart(Δt) + t, tend(Δt) + t))
end

@commutative function +(Δt::TimeInterval{<:IA.Interval}, t::Real)
    return TimeInterval(Δt.i + t)
end

function +(Δt::TimeInterval{<:Interval}, Δs::TimeInterval{<:Interval})
    return TimeInterval(minkowski_sum(Δt.i, Δs.i))
end

@commutative function +(Δt::TimeInterval{<:IA.Interval}, Δs::TimeInterval{<:Interval})
    return TimeInterval(Δt.i.dat + Δs.i)
end

function +(Δt::TimeInterval{<:IA.Interval}, Δs::TimeInterval{<:IA.Interval})
    return TimeInterval(Δt.i + Δs.i)
end

function -(Δt::TimeInterval{<:Interval}, t::Real)
    return TimeInterval(Interval(tstart(Δt) - t, tend(Δt) - t))
end

function -(Δt::TimeInterval{<:IA.Interval}, t::Real)
    return TimeInterval(Δt.i - t)
end

function isdisjoint(Δt::TimeInterval, Δs::TimeInterval)
    return tstart(Δt) > tend(Δs) || tstart(Δs) > tend(Δt)
end

function issubset(Δt::TimeInterval, Δs::TimeInterval)
    return _issubset(Δt.i, Δs.i)
end

function _issubset(x::Interval, y::Interval)
    return x ⊆ y
end

@commutative function _issubset(x::IA.Interval, y::Interval)
    return _issubset(x, y.dat)
end

function _issubset(x::IA.Interval, y::IA.Interval)
    return x ⊆ y
end

function intersect(Δt::TimeInterval, Δs::TimeInterval)
    return TimeInterval(_intersection(Δt.i, Δs.i))
end

function _intersection(x::Interval, y::Interval)
    return intersection(x, y)
end

@commutative function _intersection(x::IA.Interval, y::Interval)
    return _intersection(x, y.dat)
end

function _intersection(x::IA.Interval, y::IA.Interval)
    l, h = max(inf(x), inf(y)), min(sup(x), sup(y))
    if l > h
        return emptyinterval()
    end
    return IA.interval(l, h)
end
