"""
    post(alg::GLGM06, problem, time_horizon)
"""
function post(alg::GLGM06, problem, time_horizon)

    # ==================================
    # Initialization and discretization
    # ==================================
    max_order = alg.max_order
    δ = alg.δ
    T = time_horizon
    N = round(Int, T / δ) # number of reach sets to be evaluated at most

    # compute and unrwap discretized system
    Pdiscr = discretize(problem, δ, algorithm=alg.discretization,
                                    sih_method=alg.sih_method,
                                    exp_method=alg.exp_method,
                                    set_operations="zonotope")

    Ω0, Φ = Pdiscr.x0, Pdiscr.s.A

    # =====================
    # Flowpipe computation
    # =====================

    # preallocate output
    Rsets = Vector{ReachSet{Zonotope{Float64}}}(undef, N)

    info("Reachable States Computation...")
    @timing begin
    if inputdim(Pdiscr) == 0
        reach_homog!(Rsets, Ω0, Φ, N, δ, max_order)
    else
        U = inputset(Pdiscr) # TODO dispatch for constant input..
        reach_inhomog!(Rsets, Ω0, U, Φ, N, δ, max_order)
    end
    end # timing

    return Rsets
end
