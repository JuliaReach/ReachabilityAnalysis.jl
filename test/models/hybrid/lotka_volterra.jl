@taylorize function lotka_volterra!(du, u, p, t)
    u1u2 = u[1] * u[2]
    du[1] = 3.0 * (u[1] - u1u2)
    du[2] = u1u2 - u[2]
    return du
end

function lotka_volterra_hybrid(; ε_ext=1e-4, # threshold for the outer approximation
                                 n_int=50)   # number of directions for the inner approximation

    # generate external / internal polytopic approximations of the guard
    B = Ball2([1.0, 1.0], 0.15) # "exact"
    B_ext = overapproximate(B, ε_ext) # outer approximation
    B_int = underapproximate(B, PolarDirections(n_int)) # inner approximation
    B_int = tohrep(convert(VPolygon, B_int)) # cast to Hrep
    B_intᶜ = complement(B_int)

    # define modes
    aut = LightAutomaton(3)
    outside = @system(x' = lotka_volterra!(x), dim: 2, x ∈ B_intᶜ)
    inside = @system(x' = lotka_volterra!(x), dim: 2, x ∈ B_ext)
    outside_unconstrained = @system(x' = lotka_volterra!(x), dim: 2, x ∈ Universe(2))

    # define the transition graph
    add_transition!(aut, 1, 2, 1)
    add_transition!(aut, 2, 3, 2)
    T_out_in = @map(x -> x, dim:2, x ∈ B_ext)
    T_in_out = @map(x -> x, dim:2, x ∈ B_intᶜ)

    # initial-value problem
    H = HybridSystem(automaton=aut, modes=[outside, inside, outside_unconstrained],
                                           resetmaps=[T_out_in, T_in_out])

    # initial states
    ε = 0.008
    X0 = Hyperrectangle(low=[1.3-ε, 1.], high=[1.3+ε, 1.])

    return InitialValueProblem(H, [(X0, 1)])
end

#=
# TEST RUN:

prob = lotka_volterra_hybrid();

@time sol = solve(prob,
                  tspan=(0.0, 3.64),
                  alg=TMJets(abs_tol=1e-16, orderT=7, orderQ=1, disjointness=RA.Dummy()),
                  max_jumps=2,
                  intersect_source_invariant=false,
                  intersection_method=RA.BoxIntersection(),
                  clustering_method=RA.BoxClustering(),
                  disjointness_method=RA.BoxEnclosure());
solz = overapproximate(sol, Zonotope);

fig = plot()
plot!(fig, Ball2([1.0, 1.0], 0.15), 1e-4)
plot!(fig, solz[1], vars=(1, 2))
plot!(fig, solz[2], vars=(1, 2))
plot!(solz[3][1], vars=(1, 2))
=#
