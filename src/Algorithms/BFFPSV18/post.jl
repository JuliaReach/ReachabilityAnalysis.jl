function post(alg::BFFPSV18{N, ST}, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              time_shift::N=zero(N), kwargs...) where {N, ST}

    @unpack δ, approx_model, vars, block_indices,
            row_blocks, column_blocks = alg

    if haskey(kwargs, :NSTEPS)
        NSTEPS = kwargs[:NSTEPS]
        T = NSTEPS * δ
    else
        # get time horizon from the time span imposing that it is of the form (0, T)
        T = _get_T(tspan, check_zero=true, check_positive=true)
        NSTEPS = ceil(Int, T / δ)
    end

    # normalize system to canonical form
    # x' = Ax, x in X, or
    # x' = Ax + u, x in X, u in U
    ivp_norm = _normalize(ivp)

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, approx_model)
    Φ = state_matrix(ivp_discr)
    Ω0 = initial_state(ivp_discr)
    X = stateset(ivp_discr)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp_discr)

    # decompose the initial states into a cartesian product
    # TODO add option to do the lazy decomposition
    Xhat0 = _decompose(Ω0, column_blocks, ST)
    Φ = state_matrix(ivp_discr)
    X = stateset(ivp_discr) # invariant

    # force using sparse type for the matrix exponential
    if alg.sparse
        Φ = SparseArrays.sparse(Φ)
    end

    # variables that are actually computed
    vars = reduce(vcat, alg.row_blocks)

    # preallocate output flowpipe
    CP = CartesianProductArray{N, ST}
    F = Vector{SparseReachSet{N, CP, length(vars)}}(undef, NSTEPS)

    # option to use array views
    viewval = Val(alg.view)

    if got_homogeneous
        reach_homog_BFFPSV18!(F, Xhat0, Φ, NSTEPS, δ, X, ST,
                              vars, block_indices,
                              row_blocks, column_blocks, time_shift, viewval)

    else
        U = inputset(ivp_discr)
        @assert isa(U, LazySet) "expected input of type `<:LazySet`, but got $(typeof(U))"
        reach_inhomog_BFFPSV18!(F, Xhat0, Φ, NSTEPS, δ, X, U, ST,
                                vars, block_indices,
                                row_blocks, column_blocks, time_shift, viewval)
    end
    return Flowpipe(F)
end
