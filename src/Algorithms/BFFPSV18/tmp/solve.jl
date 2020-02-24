using SparseArrays

# x' = Ax
function solve_BFFPSV18(P::IVP{<:LCS, <:LazySet}, opts)

    # unwrap some options
    vars = opts[:vars]
    δ = opts[:δ]
    N = opts[:N]
    NUM = opts[:num_type]
    ST = opts[:set_type]
    partition = opts[:partition]
    block_indices = opts[:block_indices]
    column_blocks = opts[:column_blocks]
    row_blocks = opts[:row_blocks]

    # normalize and discretize system
    Pdiscr = discretization(P, δ)

    Ω0 = Pdiscr.x0 # bloated initial set
    # decompose into a cartesian product
    Ω0deco = _decompose(Ω0, partition, ST)
    ϕ = Pdiscr.s.A
    Xhat0 = Ω0deco.array # pointer to the CPA array

    # preallocate output flowpipe
    SRS = SparseReachSet{CartesianProductArray{NUM, ST}}
    res = Vector{SRS}(undef, N)

    # compute flowpipe, DENSE
    reach_homog!(res, ϕ, Xhat0, δ, N, vars, block_indices, row_blocks, column_blocks, NUM, ST)

    return res
end

# x' = Ax + Bu, u in U
function solve_BFFPSV18(P::IVP{<:CLCCS, <:LazySet}, opts)

    # unwrap some options
    vars = opts[:vars]
    δ = opts[:δ]
    N = opts[:N]
    NUM = opts[:num_type]
    ST = opts[:set_type]
    partition = opts[:partition]
    block_indices = opts[:block_indices]
    column_blocks = opts[:column_blocks]
    row_blocks = opts[:row_blocks]
    sp = opts[:sparse]

    # normalize and discretize system
    Pdiscr = discretization(P, δ)

    Ω0 = Pdiscr.x0 # bloated initial set
    # decompose into a cartesian product
    Ω0deco = _decompose(Ω0, partition, ST)
    ϕ = Pdiscr.s.A
    Xhat0 = Ω0deco.array # pointer to the CPA array

    # preallocate output flowpipe
    SRS = SparseReachSet{CartesianProductArray{NUM, ST}}
    res = Vector{SRS}(undef, N)

    # compute flowpipe
    U = inputset(Pdiscr.s).U # we are assuming that this input is CONSTANT and the system
    # has been normalized to be of the form x' = Ax + u, u in U

    if sp
        # SPARSE
        ϕ = sparse(ϕ) # FIXME: where is this done in Reachability.jl ?
        reach_inhomog_sparse!(res, ϕ, Xhat0, U, δ, N, vars, block_indices, row_blocks, column_blocks, NUM, ST)
    else
        # DENSE
        reach_inhomog!(res, ϕ, Xhat0, U, δ, N, vars, block_indices, row_blocks, column_blocks, NUM, ST)
    end
    return res
end
