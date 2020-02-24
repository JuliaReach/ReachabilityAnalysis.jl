export BFFPS19

# ===============================================================
# Bogomolov, Forets, Frehse, Potomkin, Schilling.
# ===============================================================

function options_BFFPS19()
    return OptionSpec[
        # general options
        OptionSpec(:discretization, "forward", domain=String,
            info="model for bloating/continuous time analysis"),
        OptionSpec(:δ, 1e-2, domain=Float64, aliases=[:sampling_time],
            domain_check=(v  ->  v > 0.), info="time step"),
        OptionSpec(:vars, Int[], domain=AbstractVector{Int}, domain_check=(
            v  ->  length(v) > 0 && all(e -> e > 0, v)),
            info="variables of interest; default: all variables"),
        OptionSpec(:guards_proj, Vector(), domain=AbstractVector,
             info="internal storage of projected guard constraints"),
        OptionSpec(:partition, [Int[]],
            domain=AbstractVector{<:AbstractVector{Int}}, domain_check=
            ispartition,
            info="block partition; a block is represented by a vector " *
                 "containing its indices"),

        # discretization options
        OptionSpec(:sih_method, "concrete", domain=String,
            info="method to compute the symmetric interval hull in discretization"),
        OptionSpec(:exp_method, "base", domain=String,
            info="method to compute the matrix exponential"),
        OptionSpec(:assume_sparse, false, domain=Bool,
            info="use an analysis for sparse discretized matrices?"),

        # reachability options
        OptionSpec(:lazy_inputs_interval, lazy_inputs_interval_always,
            domain=Union{Int, Function},
            domain_check=(v  ->  !(v isa Int) || v >= -1),
            info="length of interval in which the inputs are handled as a " *
                 "lazy set (``-1`` for 'never'); may generally also be a " *
                 "predicate over indices; the default corresponds to ``-1``"),

        # approximation options
        OptionSpec(:block_options, nothing, domain=Any,
            info="short hand to set ':block_options_init' and " *
                 "':block_options_iter'"),
        OptionSpec(:block_options_init, nothing, domain=Any,
            info="option for the approximation of the initial states " *
                 "(during decomposition)"),
        OptionSpec(:block_options_iter, nothing, domain=Any,
            info="option for the approximation of the states ``X_k``, ``k>0``"),

        # convenience options
        OptionSpec(:assume_homogeneous, false, domain=Bool,
            info="ignore dynamic inputs during the analysis?"),
    ]
end

function normalization_BFFPS19!(𝑂::TwoLayerOptions)
    # :lazy_inputs_interval option: convert integers to functions
    if haskey_specified(𝑂, :lazy_inputs_interval)
        v = 𝑂[:lazy_inputs_interval]
        if v isa Int
            if v == -1
                𝑂.specified[:lazy_inputs_interval] = lazy_inputs_interval_never
            elseif v == 0
                𝑂.specified[:lazy_inputs_interval] = lazy_inputs_interval_always
            else
                𝑂.specified[:lazy_inputs_interval] = (k -> k % v == 0)
            end
        end
    end

    # :block_options options: use convenience option for '_init' and '_iter'
    if !haskey_specified(𝑂, :block_options_init) &&
            haskey_specified(𝑂, :block_options)
        𝑂.specified[:block_options_init] = 𝑂[:block_options]
    end
    if !haskey_specified(𝑂, :block_options_iter) &&
            haskey_specified(𝑂, :block_options)
        𝑂.specified[:block_options_iter] = 𝑂[:block_options]
    end

    nothing
end

function validation_BFFPS19(𝑂)
    # block_options
    for b_options in [:block_options, :block_options_init, :block_options_iter]
        if haskey_specified(𝑂, b_options)
            bo = 𝑂[b_options]
            if bo isa Type{<:LazySet} || bo isa Type{<:AbstractDirections}
                # uniform options
            elseif bo isa Symbol
                # template directions
                option = get(Utils.template_direction_symbols, bo, nothing)
                if option == nothing
                    throw(DomainError(bo, "if the `$b_options` option " *
                        "is a Symbol, it must be one of " *
                        "$(keys(Utils.template_direction_symbols))"))
                end
                𝑂.specified[b_options] = option
            elseif bo isa AbstractVector || bo isa Dict{Int, Any}
                # mapping
            elseif bo isa Real || bo isa Pair{<:UnionAll, <:Real}
                ε = bo isa Real ? bo : bo[2]
                if ε <= 0
                    throw(DomainError(ε, "the `$b_options` option must be " *
                                         "positive"))
                end
            elseif b_options == :block_options_iter && bo == nothing
                # no overapproximation
            else
                throw(DomainError(bo == nothing ? "nothing" : bo,
                    "the `$b_options` option does not accept the given input"))
            end
        end
    end

    nothing
end

"""
    BFFPS19 <: AbstractContinuousPost

Implementation of the reachability algorithm for purely continuous linear
time-invariant systems using block decompositions by S. Bogomolov, M. Forets,
G. Frehse, A. Podelski, C. Schilling and F. Viry [1].

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Notes

The following options are available:

```julia
$(print_option_spec(options_BFFPS19()))
```

### Algorithm

We refer to [1] for technical details.

[1] [ReachabilityAnalysis analysis of linear hybrid systems via block
decomposition](https://arxiv.org/abs/1905.02458).
Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling.
"""
struct BFFPS19 <: AbstractContinuousPost
    options::TwoLayerOptions

    function BFFPS19(𝑂::Options)
        normalized_𝑂 = validate_and_wrap_options(𝑂, options_BFFPS19();
            validation=validation_BFFPS19,
            normalization=normalization_BFFPS19!)
        return new(normalized_𝑂)
    end
end

# convenience constructor from pairs of symbols
BFFPS19(𝑂::Pair{Symbol,<:Any}...) = BFFPS19(Options(Dict{Symbol,Any}(𝑂)))

# default options
BFFPS19() = BFFPS19(Options())

include("reach.jl")
include("reach_blocks.jl")

init(𝒫::BFFPS19, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝑆, copy(𝑂))

function init!(𝒫::BFFPS19, 𝑆::AbstractSystem, 𝑂::Options)
    # state dimension for (purely continuous or purely discrete systems)
    𝑂copy = copy(𝑂)
    𝑂copy[:n] = statedim(𝑆)

    # solver-specific options (adds default values for unspecified options)
    𝑂validated = validate_solver_options_and_add_default_values!(𝑂copy)

    # :vars option; default: all variables
    if haskey_specified(𝒫.options, :vars)
        𝑂validated[:vars] = 𝒫.options[:vars]
    else
        𝑂validated[:vars] = 1:𝑂validated[:n]
    end

    # :partition option: use 1D blocks
    if haskey_specified(𝒫.options, :partition)
        𝑂validated[:partition] = 𝒫.options[:partition]
    else
        𝑂validated[:partition] = [[i] for i in 1:𝑂validated[:n]]
    end

    opD = 𝒫.options[:opD]
    @assert opD isa DecomposedDiscretePost "this continuous post operator works " *
                                       "only with DecomposedDiscretePost"
    HS = 𝒫.options[:HS]
    constrained_dims = constrained_dimensions(HS)
    out_vars = opD.options[:out_vars]
    loc_id = 𝒫.options[:loc_id]
    temp_vars = unique([out_vars; constrained_dims[loc_id]])
    temp_vars = get_variables_from_relevant_blocks(𝑂validated[:partition], temp_vars)
    opD.options[:temp_vars] = temp_vars
    guards_constraints = [guard(HS, trans) for trans in out_transitions(HS, loc_id)]
    𝑂validated[:vars] = temp_vars
    𝑂validated[:guards_proj] = [project(c, temp_vars) for c in guards_constraints]
    𝑂validated[:blocks] = compute_blocks(𝑂validated[:vars], 𝑂validated[:partition])

    # :block_options_init & :block_options_iter options:
    # set default according to :partition
    if !haskey_specified(𝒫.options, :block_options_init)
        𝑂validated[:block_options_init] =
            compute_default_block_options(𝑂validated[:partition])
    end
    if !haskey_specified(𝒫.options, :block_options_iter)
        𝑂validated[:block_options_iter] =
            compute_default_block_options(𝑂validated[:partition])
    end

    if 𝑂validated[:project_reachset]
        𝑂validated[:output_function] = nothing
    else
        𝑂validated[:output_function] = 𝑂validated[:projection_matrix]
    end

    return 𝑂validated
end

"""
    post(𝒫::BFFPS19, 𝑆::AbstractSystem, 𝑂::Options)

Calculate the reachable states of the given initial value problem using `BFFPS19`.

### Input

- `𝒫` -- post operator of type `BFFPS19`
- `𝑆` -- sytem, initial value problem for a continuous ODE
- `𝑂` -- algorithm-specific options
"""
function post(𝒫::BFFPS19, 𝑆::AbstractSystem, 𝑂_input::Options)
    𝑂 = TwoLayerOptions(merge(𝑂_input, 𝒫.options.specified), 𝒫.options.defaults)

    @assert 𝑂[:mode] == "reach" "the mode $(𝑂[:mode]) is not supported"

    info("Reachable States Computation...")
    @timing begin
        Rsets = reach_mixed(𝑆, 𝑂)

        info("- Total")
    end

    # Projection
    if 𝑂[:project_reachset]
        info("Projection...")
        RsetsProj = @timing project(Rsets, 𝑂)
    else
        RsetsProj = Rsets
    end

    return ReachSolution(RsetsProj, 𝑂_input)
end

function get_variables_from_relevant_blocks(partition, vars)
    result = Vector{Int}()
    for ur in partition
        if any(v ∈ ur for v in vars)
            append!(result, collect(ur))
        end
    end
    return result
end
