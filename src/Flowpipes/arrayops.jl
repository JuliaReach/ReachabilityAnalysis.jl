# ========================================================================
# Specialized functions on SingleEntryVector
# TODO refactor to LazySets.jl
# ========================================================================

# difference between single entry vectors (SEV); it is type unstable because
# it may return either a SEV or a regular Vector depending on the indices
# but it is much faster than the fallback if the indices match
function minus(e1::SingleEntryVector{N}, e2::SingleEntryVector{N}) where {N}
    e1.n == e2.n ||
        throw(DimensionMismatch("dimensions must match, but they are $(length(e1)) and $(length(e2)) respectively"))

    if e1.i == e2.i
        return SingleEntryVector(e1.i, e1.n, e1.v - e2.v)
    else
        out = zeros(N, e1.n)
        @inbounds begin
            out[e1.i] = e1.v
            out[e2.i] = -e2.v
        end
        return out
    end
end

# addition between single entry vectors (SEV); it is type unstable because
# it may return either a SEV or a regular Vector depending on the indices
# but it is much faster than the fallback if the indices match
function plus(e1::SingleEntryVector{N}, e2::SingleEntryVector{N}) where {N}
    e1.n == e2.n ||
        throw(DimensionMismatch("dimensions must match, but they are $(length(e1)) and $(length(e2)) respectively"))

    if e1.i == e2.i
        return SingleEntryVector(e1.i, e1.n, e1.v + e2.v)
    else
        out = zeros(N, e1.n)
        @inbounds begin
            out[e1.i] = e1.v
            out[e2.i] = e2.v
        end
        return out
    end
end

# norm of the difference of two SEV ||x - y|| in the (vector p norm)
function normdiff(e1::SingleEntryVector{N}, e2::SingleEntryVector{N}, p::Real=2) where {N}
    e1.n == e2.n ||
        throw(DimensionMismatch("dimensions must match, but they are $(length(e1)) and $(length(e2)) respectively"))

    if e1.i == e2.i
        δ = e1.v - e2.v
        return abs(δ)

    else
        a = abs(e1.v)
        b = abs(e2.v)
        if isinf(p)
            return max(a, b)
        else
            s = a^p + b^p
            return s^(1 / p)
        end
    end
end

# fallback
function normdiff(x::AbstractVector{N}, y::AbstractVector{N}, p::Real=2) where {N}
    return norm(x - y, p)
end
