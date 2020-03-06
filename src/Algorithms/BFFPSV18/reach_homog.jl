# ====================
# x' = Ax
# (no input)
# ====================

# Set representation: Interval
# Matrix operations: Dense
# Invariant: No
function reach_homog!(F, Xhat0, Φ, NSTEPS, δ, ::Universe, vars, block_indices,
                      row_blocks, column_blocks, N, ST::Type{<:Interval})

    # store first element
    R0 = CartesianProductArray(Xhat0[block_indices])
    Δt = zero(N) .. δ
    res[1] = SparseReachSet(R0, Δt, vars) # TODO : use view?

    # cache matrix
    ϕpowerk = copy(ϕ)
    ϕpowerk_cache = similar(ϕ)

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))
    zero_int = Interval(zero(N), zero(N))

    @inbounds for k in 2:N
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            Xhatk[i] = zero_int
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                # concrete interval-interval Minkowski sum
                aux = ϕpowerk[bi[1], bj[1]] * Xhat0[j].dat # TODO : concrete linear_map from LazySets ?
                Xhatk[i] = minkowski_sum(Xhatk[i], Interval(aux))
            end
        end
        Δt += δ
        Rk = CartesianProductArray(copy(Xhatk))
        res[k] = SparseReachSet(Rk, Δt, vars)

        mul!(ϕpowerk_cache, ϕpowerk, ϕ)
        copyto!(ϕpowerk, ϕpowerk_cache)
    end
    return res
end

# generic setrep, no input
# assuming that Xhat0 Xhat0 is already decomposed
function _reach_homog!(res, ϕ, Xhat0, δ, N, vars, block_indices, row_blocks, column_blocks, NUM, ST)

    # store first reachset
    ti, tf = 0.0, δ
    dt = 0 .. δ
    R0 = ReachSet(CartesianProductArray(Xhat0[block_indices]), ti, tf)
    res[1] = SparseReachSet(R0, vars)

    # cache matrix
    ϕpowerk = copy(ϕ)
    ϕpowerk_cache = similar(ϕ)

    # preallocate buffer for each row-block
    buffer = Vector{LazySets.LinearMap{NUM, ST, NUM, Matrix{NUM}}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    @inbounds for k in 2:N
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                buffer[j] = ϕpowerk[bi, bj] * Xhat0[j]
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer), ST)
        end
        ti = tf
        tf += δ
        Rk = ReachSet(CartesianProductArray(copy(Xhatk)), ti, tf)
        res[k] = SparseReachSet(Rk, vars)

        mul!(ϕpowerk_cache, ϕpowerk, ϕ)
        copyto!(ϕpowerk, ϕpowerk_cache)
    end
    return res
end

# generic setrep, no input, sparse
function _reach_homog_sparse!(res, ϕ, Xhat0, δ, N, vars, block_indices, row_blocks, column_blocks, NUM, ST)

    # store first element
    ti, tf = 0.0, δ
    R0 = ReachSet(CartesianProductArray(Xhat0[block_indices]), ti, tf)
    res[1] = SparseReachSet(R0, vars)

    # cache matrix
    ϕpowerk = copy(ϕ)
    ϕpowerk_cache = similar(ϕ)

    # preallocate buffer for each row-block
    buffer = Vector{LazySets.LinearMap{NUM, ST, NUM, Matrix{NUM}}}(undef, length(column_blocks))

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    @inbounds for k in 2:N
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                ϕpowerk_bi_bj = ϕpowerk
                if !iszero(ϕpowerk_bi_bj)
                    buffer[j] = ϕpowerk[bi, bj] * Xhat0[j]
                end
            end
            Xhatk[i] = overapproximate(MinkowskiSumArray(buffer), ST)
        end
        ti = tf
        tf += δ
        Rk = ReachSet(CartesianProductArray(copy(Xhatk)), ti, tf)
        res[k] = SparseReachSet(Rk, vars)

        mul!(ϕpowerk_cache, ϕpowerk, ϕ)
        copyto!(ϕpowerk, ϕpowerk_cache)
    end
    return res
end
