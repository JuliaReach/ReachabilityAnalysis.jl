function options_ASB07_decomposed()
    𝑂spec = Vector{OptionSpec}()

    push!(𝑂spec, OptionSpec(:δ, 1e-2, domain=Float64, aliases=[:sampling_time],
        domain_check=(v  ->  v > 0.), info="time step"))

    push!(𝑂spec, OptionSpec(:max_order, 10, domain=Int,
        info="maximum allowed order of zonotopes"))

    push!(𝑂spec, OptionSpec(:order_discretization, 2, domain=Int,
        info="order of Taylor approximation in discretization"))

    push!(𝑂spec, OptionSpec(:set_operations_discretization, "zonotope",
        domain=String, domain_check=(v  ->  v ∈ ["zonotope", "lazy"]),
        info="type of set operations applied during discretization"))

    push!(𝑂spec, OptionSpec(:block_options_init, Hyperrectangle, domain=Any,
        info="option for the decomposition of the initial states"))

    push!(𝑂spec, OptionSpec(:partition, [Int[]],
        domain=AbstractVector{<:AbstractVector{Int}}, domain_check=ispartition,
        info="block partition; a block is represented by a vector containing " *
             "its indices"))

    push!(𝑂spec, OptionSpec(:vars, Int[], domain=AbstractVector{Int},
        domain_check=(v  ->  length(v) > 0 && all(e -> e > 0, v)),
        info="variables of interest; default: all variables"))

    return 𝑂spec
end
