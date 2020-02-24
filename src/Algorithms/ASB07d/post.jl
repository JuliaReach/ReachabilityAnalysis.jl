function post(𝒜::ADB07d,
              𝑃::InitialValueProblem{<:AbstractContinuousSystem},
              𝑂::Options)
    # =================================
    # Initialization and discretization
    # =================================

    𝑂 = merge(𝒜.options.defaults, 𝑂, 𝒜.options.specified)
    δ, T = 𝑂[:δ], 𝑂[:T]
    N = round(Int, T / δ)
    n = 𝑂[:n]
    partition = 𝑂[:partition]
    blocks = 𝑂[:blocks]
    max_order = 𝑂[:max_order]

    # compute and unrwap discretized system
    𝑃_discrete = discretize(𝑃, δ; algorithm="interval_matrix",
                            order=𝑂[:order_discretization],
                            set_operations=𝑂[:set_operations_discretization])
    Ω0, Φ = 𝑃_discrete.x0, 𝑃_discrete.s.A

    # decompose initial states
    # TODO this should use the concrete linear_map projection (see LazySets#1726)
    Ωhat0 = array(decompose(Ω0, partition, 𝑂[:block_options_init]))

    # ====================
    # Flowpipe computation
    # ====================

    # preallocate output
    Rsets = Vector{ReachSet{CartesianProductArray{Float64, LazySet{Float64}}}}(undef, N)

    info("Reachable States Computation...")
    @timing begin
    if inputdim(𝑃_discrete) == 0
        U = nothing
    else
        U = inputset(𝑃_discrete)
    end
    reach_ASB07_decomposed!(Rsets, Ωhat0, U, Φ, N, δ, max_order, n, partition,
                            blocks)
    end # timing

    Rsol = ReachSolution(Rsets, 𝑂)

    # ==========
    # Projection
    # ==========

    if 𝑂[:project_reachset]
        info("Projection...")
        Rsol = @timing project(Rsol)
    end

    return Rsol
end
