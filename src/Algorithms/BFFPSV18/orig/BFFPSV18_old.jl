export BFFPSV18

# TODO: DOCS:
# see branch
# https://github.com/JuliaReach/Reachability.jl/pull/499
# see branch mforets/decompose_options

# ===============================================================
# Bogomolov, Forets, Frehse, Podelski, Schilling, Viry. HSCC 2018
# ===============================================================

# dummy functions for option :lazy_inputs_interval
lazy_inputs_interval_always = (k -> true)
lazy_inputs_interval_never = (k -> false)

function ispartition(partition::AbstractVector{<:AbstractVector{Int}})
    current = 1
    for block in partition
        for i in block
            if i != current
                return false
            end
            current += 1
        end
    end
    return true
end

function options_BFFPSV18()
    return OptionSpec[
        # general options
        OptionSpec(:discretization, "forward", domain=String,
            info="model for bloating/continuous time analysis"),
        OptionSpec(:algorithm, "explicit", domain=String, domain_check=(
            v  ->  v in ["explicit", "wrap"]), info="algorithm backend"),
        OptionSpec(:δ, 1e-2, domain=Float64, aliases=[:sampling_time],
            domain_check=(v  ->  v > 0.), info="time step"),
        OptionSpec(:vars, Int[], domain=AbstractVector{Int}, domain_check=(
            v  ->  length(v) > 0 && all(e -> e > 0, v)),
            info="variables of interest; default: all variables"),
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
        OptionSpec(:eager_checking, true, domain=Bool,
            info="terminate as soon as property violation was detected?"),
    ]
end

function normalization_BFFPSV18!(𝑂::TwoLayerOptions)
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

function validation_BFFPSV18(𝑂)
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
    BFFPSV18 <: AbstractContinuousPost

Implementation of the reachability algorithm for purely continuous linear
time-invariant systems using block decompositons by S. Bogomolov, M. Forets,
G. Frehse, A. Podelski, C. Schilling and F. Viry [1].

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Notes

The following options are available:

```julia
$(print_option_spec(options_BFFPSV18()))
```

### Algorithm

We refer to [1] for technical details.

[1] [Reach Set Approximation through Decomposition with Low-dimensional Sets
and High-dimensional Matrices](https://dl.acm.org/citation.cfm?id=3178128).
S. Bogomolov, M. Forets, G. Frehse, A. Podelski, C. Schilling, F. Viry.
HSCC '18 Proceedings of the 21st International Conference on Hybrid Systems:
Computation and Control (part of CPS Week).
"""
struct BFFPSV18 <: AbstractContinuousPost
    options::TwoLayerOptions

    function BFFPSV18(𝑂::Options)
        normalized_𝑂 = validate_and_wrap_options(𝑂, options_BFFPSV18();
            validation=validation_BFFPSV18,
            normalization=normalization_BFFPSV18!)
        return new(normalized_𝑂)
    end
end

# convenience constructor from pairs of symbols
BFFPSV18(𝑂::Pair{Symbol,<:Any}...) = BFFPSV18(Options(Dict{Symbol,Any}(𝑂)))

# default options
BFFPSV18() = BFFPSV18(Options())

init(𝒫::BFFPSV18, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝑆, copy(𝑂))

function init!(𝒫::BFFPSV18, 𝑆::AbstractSystem, 𝑂::Options)
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

    # :blocks option (internal only)
    # list of all interesting block indices in the partition
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
    post(𝒫::BFFPSV18, 𝑆::AbstractSystem, 𝑂::Options)

Calculate the reachable states of the given initial value problem using `BFFPSV18`.

### Input

- `𝒫` -- post operator of type `BFFPSV18`
- `𝑆` -- sytem, initial value problem for a continuous ODE
- `𝑂` -- algorithm-specific options
"""
function post(𝒫::BFFPSV18, 𝑆::AbstractSystem, 𝑂_input::Options)
    𝑂 = TwoLayerOptions(merge(𝑂_input, 𝒫.options.specified), 𝒫.options.defaults)

    if 𝑂[:mode] == "reach"
        info("Reachable States Computation...")
        @timing begin
            Rsets = reach(𝑆, 𝑂)
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

    elseif 𝑂[:mode] == "check"
        info("invariants are currently not supported in 'check' mode")

        # Input -> Output variable mapping in property
        property = inout_map_property(𝑂[:property], 𝑂[:partition], 𝑂[:blocks], 𝑂[:n])

        # =================
        # Property checking
        # =================
        info("Property Checking...")
        @timing begin
            answer = check_property(𝑆, property, 𝑂)
            info("- Total")
        end

        if answer == 0
            info("The property is satisfied!")
            return CheckSolution(true, -1, 𝑂_input)
        else
            info("The property may be violated at index $answer," *
                " (time point $(answer * 𝑂[:δ]))!")
            return CheckSolution(false, answer, 𝑂_input)
        end
    else
        error("unsupported mode $(𝑂[:mode])")
    end # mode
end

function compute_blocks(vars, partition)
    blocks = Vector{Int}()
    sizehint!(blocks, length(vars))
    next = 0
    var_idx = 1
    for (i, block) in enumerate(partition)
        next += length(block)
        if vars[var_idx] <= next
            push!(blocks, i)
            var_idx += 1
            while var_idx <= length(vars) && vars[var_idx] <= next
                var_idx += 1
            end
            if var_idx > length(vars)
                break
            end
        end
    end
    @assert var_idx == length(vars) + 1
    sizehint!(blocks, length(blocks))
    return blocks
end

function compute_default_block_options(partition)
    # use Interval for 1D blocks and Hyperrectangle otherwise
    block_options = Vector{Type{<:LazySet}}(undef, length(partition))
    for (i, block) in enumerate(partition)
        block_options[i] = length(block) == 1 ? Interval : Hyperrectangle
    end
    return block_options
end
