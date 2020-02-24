function options_ASB07()
    𝑂spec = Vector{OptionSpec}()

    push!(𝑂spec, OptionSpec(:δ, 1e-2, domain=Float64, aliases=[:sampling_time],
                            domain_check=(v  ->  v > 0.), info="time step"))

    push!(𝑂spec, OptionSpec(:max_order, 10, domain=Int,
                            info="maximum allowed order of zonotopes"))

    push!(𝑂spec, OptionSpec(:order_discretization, 2, domain=Int,
                            info="order of Taylor approximation in " *
                                 "discretization"))

    push!(𝑂spec, OptionSpec(:set_operations_discretization, "zonotope",
                            domain=String,
                            domain_check=(v  ->  v ∈ ["zonotope", "lazy"]),
                            info="type of set operations applied during " *
                            "discretization"))

    return 𝑂spec
end
