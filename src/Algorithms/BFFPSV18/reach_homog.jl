# =======================
# x' = Ax
# =======================

# Set representation: Interval
# Matrix operations: Dense
# Invariant: No
function reach_homog_BFFPSV18!(F, Xhat0, Φ::AbstractMatrix{NM}, NSTEPS, δ, X::Universe,
                               setrep::Type{<:Interval{N}}, vars, block_indices,
                               row_blocks, column_blocks, time_shift::N) where {NM, N}

    # initial reach set
    Δt = (zero(N) .. δ) + time_shift
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0.array[collect(block_indices)]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # buffer for overapproximated Minkowski sum of each row-block
    Xhatk = Vector{setrep}(undef, length(row_blocks))
    zero_int = Interval(zero(N), zero(N))

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            Xhatk[i] = zero_int
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                # concrete interval-interval Minkowski sum
                aux = Φpowerk[bi[1], bj[1]] * Xhat0.array[j].dat
                Xhatk[i] = Interval(Xhatk[i].dat + aux)
            end
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)
    end
    return F
end

# Set representation: Generic
# Matrix operations: Dense
# Invariant: No
function reach_homog_BFFPSV18!(F, Xhat0, Φ::AbstractMatrix{NM}, NSTEPS, δ, X::Universe,
                               setrep, vars, block_indices, row_blocks, column_blocks, time_shift::N) where {NM, N}

    # initial reach-set
    Δt = (zero(N) .. δ) + time_shift
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    buffer = Vector{LazySets.LinearMap{N, ST, N, Matrix{N}}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                buffer[j] = Φpowerk[bi, bj] * Xhat0[j]
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer), ST)
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)
    end
    return F
end

# Set representation: Generic
# Matrix operations: Sparse
# Invariant: No
function reach_homog_BFFPSV18!(F, Xhat0, Φ::SparseMatrixCSC{NM, IM}, NSTEPS, δ::N, X::Universe,
                               setrep, vars, block_indices,
                               row_blocks, column_blocks, time_shift::N) where {NM, IM, N}

    # store first element
    Δt = (zero(N) .. δ) + time_shift
    @inbounds F[1] = SparseReachSet(CartesianProductArray(Xhat0[block_indices]), Δt, vars)

    # cache matrix
    Φpowerk = copy(Φ)
    Φpowerk_cache = similar(Φ)

    # preallocate buffer for each row-block
    buffer = Vector{LazySets.LinearMap{N, ST, NM, SparseMatrixCSC{NM, IM}}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    @inbounds for k in 2:NSTEPS
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                Φpowerk_bi_bj = Φpowerk
                if !iszero(Φpowerk_bi_bj)
                    buffer[j] = Φpowerk[bi, bj] * Xhat0[j]
                end
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer), ST)
        end
        Δt += δ
        Xk = CartesianProductArray(copy(Xhatk))
        F[k] = SparseReachSet(Xk, Δt, vars)

        mul!(Φpowerk_cache, Φpowerk, Φ)
        copyto!(Φpowerk, Φpowerk_cache)
    end
    return F
end
