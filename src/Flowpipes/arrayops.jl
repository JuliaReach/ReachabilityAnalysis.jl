# ======================================================================
# Logarithmic norms of matrices (a.k.a matrix measures)
#
# Reference:
#
# C. Desoer and M. Vidyasagar, Feedback Systems: Input-Output Properties.
# Philadelphia, PA: Society for Industrial and Applied Mathematics,
# Jan. 2009, ser. Classics in Applied Mathematics.
# ======================================================================

# logarithmic norm (also known as matrix measure) for commonly used p-norms
function logarithmic_norm(A::AbstractMatrix, p::Real=Inf)
    if p == Inf
        return logarithmic_norm_inf(A)
    elseif p == 1
        return logarithmic_norm_1(A)
    elseif p == 2
        return logarithmic_norm_2(A)
    else
        throw(ArgumentError("logarithmic norm only implemented for p = 1, 2, or Inf, got p = $p"))
    end
end

# max_j a_jj + ∑_{i ≠ j} |a_ij|
function logarithmic_norm_1(A::AbstractMatrix{N}) where {N}
    out = -Inf
    @inbounds for j in 1:size(A, 2)
        α = A[j, j]
        for i in 1:size(A, 1)
            if i ≠ j
                α += abs(A[i, j])
            end
        end
        if α > out
            out = α
        end
    end
    return out
end

# max_i a_ii + ∑_{j ≠ i} |a_ij|
function logarithmic_norm_inf(A::AbstractMatrix{N}) where {N}
    out = -Inf
    @inbounds for i in 1:size(A, 1)
        α = A[i, i]
        for j in 1:size(A, 2)
            if i ≠ j
                α += abs(A[i, j])
            end
        end
        if α > out
            out = α
        end
    end
    return out
end

# max_j  1/2 * λⱼ(A + A^T)
function logarithmic_norm_2(A::AbstractMatrix)
    B = A + A'
    λ = eigvals(B)
    return maximum(λ) / 2
end

# ========================================================================
# Specialized functions on LazySets.Arrays.SingleEntryVector
# TODO refactor to LazySets.jl
# ========================================================================

# difference between single entry vectors (SEV); it is type unstable because
# it may return either a SEV or a regular Vector depending on the indices
# but it is much faster than the fallback if the indices match
function minus(e1::SingleEntryVector{N}, e2::SingleEntryVector{N}) where {N}
    e1.n == e2.n || throw(DimensionMismatch("dimensions must match, but they are $(length(e1)) and $(length(e2)) respectively"))

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
    e1.n == e2.n || throw(DimensionMismatch("dimensions must match, but they are $(length(e1)) and $(length(e2)) respectively"))

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
    e1.n == e2.n || throw(DimensionMismatch("dimensions must match, but they are $(length(e1)) and $(length(e2)) respectively"))

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
            return s^(1/p)
        end
    end
end

# fallback
function normdiff(x::AbstractVector{N}, y::AbstractVector{N}, p::Real=2) where {N}
    norm(x - y, p)
end
