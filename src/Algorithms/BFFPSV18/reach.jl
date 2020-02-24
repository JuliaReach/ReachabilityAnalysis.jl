# ====================
# x' = Ax
# ====================

# generic setrep, no input
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

# INTERVAL, no input
function reach_homog!(res, ϕ, Xhat0, δ, N, vars, block_indices, row_blocks, column_blocks, NUM, ST::Type{<:Interval})

    # store first element
    ti, tf = 0.0, δ
    R0 = ReachSet(CartesianProductArray(Xhat0[block_indices]), ti, tf)
    res[1] = SparseReachSet(R0, vars)

    # cache matrix
    ϕpowerk = copy(ϕ)
    ϕpowerk_cache = similar(ϕ)

    # preallocate overapproximated Minkowski sum for each row-block
    Xhatk = Vector{ST}(undef, length(row_blocks))

    @inbounds for k in 2:N
        for (i, bi) in enumerate(row_blocks) # loop over row-blocks of interest
            Xhatk[i] = Interval(zero(NUM), zero(NUM))
            for (j, bj) in enumerate(column_blocks) # loop over all column-blocks
                Xhatk[i] += Interval(ϕpowerk[bi[1], bj[1]] * Xhat0[j].dat)
            end
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
