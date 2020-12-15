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
