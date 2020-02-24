export DecomposedDiscretePost

"""
    DecomposedDiscretePost <: AbstractDiscretePost

Textbook implementation of a discrete post operator, but with decomposed intersections.

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Algorithm

The algorithm is based on [ReachabilityAnalysis analysis of linear hybrid systems via block
decomposition](https://arxiv.org/pdf/1905.02458.pdf).
This discrete-post operator can only be used in combination with the continuous-post operator [`BFFPS19`](@ref).
"""
struct DecomposedDiscretePost <: AbstractDiscretePost
    options::Options

    function DecomposedDiscretePost(𝑂::Options)
        𝑂copy = copy(𝑂)
        # TODO: Check why it takes always default value for convex_hull
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:overapproximation], Hyperrectangle)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:out_vars], Vector{Int}())

        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
DecomposedDiscretePost(𝑂::Pair{Symbol,<:Any}...) = DecomposedDiscretePost(Options(Dict{Symbol,Any}(𝑂)))

# default options for the DecomposedDiscretePost discrete post operator
DecomposedDiscretePost() = DecomposedDiscretePost(Options())

init(𝒫::DecomposedDiscretePost, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

# TODO: use 𝑂 only?
function init!(𝒫::DecomposedDiscretePost, 𝒮::AbstractSystem, 𝑂::Options)
    𝑂[:n] = statedim(𝒮, 1)

    # solver-specific options (adds default values for unspecified options)
    𝑂out = validate_solver_options_and_add_default_values!(𝑂)

    return 𝑂out
end

function tube⋂inv!(𝒫::DecomposedDiscretePost,
                   reach_tube::Vector{<:SparseReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  ) where {N}

    # counts the number of sets R⋂I added to Rsets
    count = 0
    @inbounds for reach_set in reach_tube
        # intersection with invariant is computed inside BFFPS19 CPost operator
        push!(Rsets,
              substitute(reach_set,
                         time_start=time_start(reach_set) + start_interval[1],
                         time_end=time_end(reach_set) + start_interval[2]))
        count = count + 1
    end

    return count
end

function post(𝒫::DecomposedDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, <:AbstractReachSet{LazySet{N}}, Int}},
              passed_list,
              source_loc_id,
              tube⋂inv,
              count_Rsets,
              jumps,
              options
             ) where {N}
    jumps += 1
    n = options[:n]
    oa = 𝒫.options[:overapproximation]
    n_lowdim = length(𝒫.options[:temp_vars])
    source_invariant = HS.modes[source_loc_id].X
    inv_isa_Hrep, inv_isa_H_polytope = get_Hrep_info(source_invariant)
    for trans in out_transitions(HS, source_loc_id)
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        target_loc = HS.modes[target(HS, trans)]
        target_invariant = target_loc.X
        constrained_map = resetmap(HS, trans)
        guard = stateset(constrained_map)
        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}}}()
        sizehint!(post_jump, count_Rsets)
        for reach_set in tube⋂inv[length(tube⋂inv) - count_Rsets + 1 : end]
            if (dim(set(reach_set)) == n_lowdim && n_lowdim < n)
                continue
            end
            # check intersection with guard
            R⋂G = Intersection(set(reach_set), guard)
            if isempty(R⋂G)
                continue
            end
            R⋂G = overapproximate(R⋂G, CartesianProductArray, oa)

            # apply assignment
            A⌜R⋂G⌟ = apply_assignment(𝒫, constrained_map, R⋂G)
            A⌜R⋂G⌟ = overapproximate(A⌜R⋂G⌟, CartesianProductArray, oa)

            # intersect with target invariant
            A⌜R⋂G⌟⋂I = Intersection(target_invariant, A⌜R⋂G⌟)
            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end

            A⌜R⋂G⌟⋂I = overapproximate(A⌜R⋂G⌟⋂I, CartesianProductArray, oa)

            # store result
            push!(post_jump, ReachSet{LazySet{N}}(A⌜R⋂G⌟⋂I,
                                                  time_start(reach_set),
                                                  time_end(reach_set)))
        end

        postprocess(𝒫, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps)
    end
end
