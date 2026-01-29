# general method given a reachability algorithm for the linear system
function reach_CARLIN_alg(X0, F1, F2; alg, resets, N, T, Δt0, bloat, compress)
    if resets == 0
        reach_CARLIN(X0, F1, F2; alg, N, T, Δt=Δt0, bloat, compress)
    else
        reach_CARLIN_resets(X0, F1, F2, resets; alg, N, T, Δt=Δt0, bloat, compress)
    end
end

function reach_CARLIN(X0, F1, F2; alg, N, T, Δt, bloat, compress, A=nothing)
    error_bound_func = error_bound_specabs

    # lift initial states
    if compress
        ŷ0 = lift_vector(X0, N) # see CarlemanLinearization/linearization.jl
    else
        ŷ0 = box_approximation(kron_pow_stack(X0, N))
    end

    # solve continuous ODE
    if isnothing(A)
        A = build_matrix(F1, F2, N; compress)
    end
    prob = @ivp(ŷ' = A * ŷ, ŷ(0) ∈ ŷ0)
    sol = solve(prob; T=T, alg=alg, Δt0=Δt)

    # projection onto the first n variables
    n = dim(X0)
    πsol_1n = _project(sol, 1:n)

    if !bloat
        return πsol_1n
    end

    # compute errors
    errfunc = error_bound_func(X0, Matrix(F1), Matrix(F2); N=N)

    # evaluate error bounds for each reach-set in the solution
    E = [errfunc.(tspan(R)) for R in sol]

    # if the interval is always > 0 then we can just take max(Ei)

    # symmetrize intervals
    E_rad = [symmetric_interval_hull(Interval(ei)) for ei in E]
    E_ball = [BallInf(zeros(n), high(ei, 1)) for ei in E_rad]

    # sum the solution with the error
    fp_bloated = Flowpipe([ReachSet(set(Ri) ⊕ Ei, tspan(Ri)) for (Ri, Ei) in zip(πsol_1n, E_ball)])

    return fp_bloated
end

function _compute_resets(resets::Int, T)
    return mince(0 .. T, resets + 1)
end

function _compute_resets(resets::Vector{Float64}, T)
    # assumes initial time is 0
    aux = vcat(0.0, resets, T)
    return [interval(aux[i], aux[i + 1]) for i in 1:(length(aux) - 1)]
end

function reach_CARLIN_resets(X0, F1, F2, resets; alg, N, T, Δt, bloat, compress)

    # build state matrix (remains unchanged upon resets)
    A = build_matrix(F1, F2, N; compress)

    # time intervals to compute
    time_intervals = _compute_resets(resets, T)

    # compute until first chunk
    T1 = sup(first(time_intervals))
    sol_1 = reach_CARLIN(X0, F1, F2; alg, N, T=T1, Δt=interval(0) + Δt, bloat, compress, A=A)

    # preallocate output flowpipe
    fp_1 = flowpipe(sol_1)
    out = Vector{typeof(fp_1)}()
    push!(out, fp_1)

    # approximate final reach-set
    Rlast = sol_1[end]
    X0 = box_approximation(set(Rlast))

    # compute remaining chunks
    for i in 2:length(time_intervals)
        T0 = T1
        Ti = sup(time_intervals[i])
        sol_i = reach_CARLIN(X0, F1, F2; alg, N, T=Ti - T0, Δt=interval(T0) + Δt, bloat, compress,
                             A=A)
        push!(out, flowpipe(sol_i))
        X0 = box_approximation(set(sol_i[end]))
        T1 = Ti
    end

    return MixedFlowpipe(out)
end
