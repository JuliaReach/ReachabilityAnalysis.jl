function post(𝒜::ASB07,
              𝑃::InitialValueProblem{<:AbstractContinuousSystem},
              𝑂::Options)
    # =================================
    # Initialization and discretization
    # =================================

    𝑂 = merge(𝒜.options.defaults, 𝑂, 𝒜.options.specified)
    δ, T = 𝑂[:δ], 𝑂[:T]
    N = round(Int, T / δ)

    # compute and unrwap discretized system
    𝑃_discrete = discretize(𝑃, δ; algorithm="interval_matrix",
                            order=𝑂[:order_discretization],
                            set_operations=𝑂[:set_operations_discretization])
    Ω0, Φ = 𝑃_discrete.x0, 𝑃_discrete.s.A

    # ====================
    # Flowpipe computation
    # ====================

    # preallocate output
    T = 𝑂[:set_operations_discretization] == "zonotope" ? Zonotope : LazySet
    Rsets = Vector{ReachSet{T{Float64}}}(undef, N)
    # Flowpipe(args.ST, N)

    max_order = 𝑂[:max_order]

    info("Reachable States Computation...")
    @timing begin
    if inputdim(𝑃_discrete) == 0
        U = nothing
    else
        U = inputset(𝑃_discrete)
    end
    reach_ASB07!(Rsets, Ω0, U, Φ, N, δ, max_order)
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
