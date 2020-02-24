using LazySets: CachedMinkowskiSumArray

"""
    check_property(S, property, options)

Interface to property checking algorithms for an LTI system.

### Input

- `S`        -- LTI system, discrete or continuous
- `property` -- property
- `options`  -- additional options

### Notes

A dictionary with available algorithms is available via
`available_algorithms_check`.
"""
function check_property(S::IVP{<:AbstractDiscreteSystem},
                        property::Property,
                        options::TwoLayerOptions
                       )::Int
    # optional matrix conversion
    S = matrix_conversion(S, options)

    # list containing the arguments passed to any reachability function
    args = []

    # coefficients matrix
    A = S.s.A
    push!(args, A)

    # determine analysis mode (sparse/dense) for lazy_expm mode
    if A isa SparseMatrixExp
        push!(args, Val(options[:assume_sparse]))
    end

    n = statedim(S)
    blocks = options[:blocks]
    partition = convert_partition(options[:partition])
    T = options[:T]
    N = (T == Inf) ? nothing : ceil(Int, T / options[:δ])

    # Cartesian decomposition of the initial set
    if length(partition) == 1 && length(partition[1]) == n &&
            options[:block_options_init] == LinearMap
        info("- Skipping decomposition of X0")
        Xhat0 = [S.x0]
    else
        info("- Decomposing X0")
        @timing begin
            Xhat0 = array(decompose(S.x0, options[:partition],
                                    options[:block_options_init]))
        end
    end

    # shortcut if only the initial set is required
    if N == 1
        if length(blocks) == 1
            Xhat0_mod = Xhat0[blocks[1]]
        else
            Xhat0_mod = CartesianProductArray(Xhat0)
        end
        return check(property, Xhat0_mod) ? 0 : 1
    end
    push!(args, Xhat0)

    # inputs
    if !options[:assume_homogeneous] && inputdim(S) > 0
        U = inputset(S)
    else
        U = nothing
    end
    push!(args, U)

    # overapproximation function for states
    block_options_iter = options[:block_options_iter]
    if block_options_iter isa AbstractVector ||
            block_options_iter isa Dict{Int, Any}
        # individual overapproximation options per block
        overapproximate_fun = (i, X) -> overapproximate(X, block_options_iter[i])
    else
        # uniform overapproximation options for each block
        overapproximate_fun = (i, X) -> overapproximate(X, block_options_iter)
    end

    # overapproximate function for inputs
    lazy_inputs_interval = options[:lazy_inputs_interval]
    if lazy_inputs_interval == lazy_inputs_interval_always
        overapproximate_inputs_fun = (k, i, x) -> overapproximate_fun(i, x)
    else
        # first set in a series
        function _f(k, i, x::LinearMap{NUM}) where {NUM}
            @assert k == 1 "a LinearMap is only expected in the first iteration"
            return CachedMinkowskiSumArray(LazySet{NUM}[x])
        end
        # further sets of the series
        function _f(k, i, x::MinkowskiSum{NUM, <:CachedMinkowskiSumArray}) where NUM
            if has_constant_directions(block_options_iter, i)
                # forget sets if we do not use epsilon-close approximation
                forget_sets!(x.X)
            end
            push!(array(x.X), x.Y)
            if lazy_inputs_interval(k)
                # overapproximate lazy set
                y = overapproximate_fun(i, x.X)
                return CachedMinkowskiSumArray(LazySet{NUM}[y])
            end
            return x.X
        end
        function _f(k, i, x)
            # other set types
            if lazy_inputs_interval(k)
                # overapproximate lazy set
                return overapproximate_fun(i, x.X)
            end
            return x
        end
        overapproximate_inputs_fun = _f
    end
    push!(args, overapproximate_inputs_fun)

    # ambient dimension
    push!(args, n)

    # number of computed sets
    push!(args, N)

    # add mode-specific block(s) argument
    algorithm = options[:algorithm]
    if algorithm == "explicit"
        push!(args, blocks)
        push!(args, partition)
        algorithm_backend = "explicit_blocks"
    else
        error("Unsupported algorithm: ", algorithm)
    end

    # add eager/lazy checking option
    push!(args, options[:eager_checking])

    # add property
    push!(args, property)

    info("- Computing successors")

    # progress meter
    progress_meter = (N != nothing) ?
        Progress(N, 1, "Computing successors ") :
        nothing
    push!(args, progress_meter)

    # call the adequate function with the given arguments list
    answer =
        @timing available_algorithms_check[algorithm_backend]["func"](args...)

    # return the result
    return answer
end

function check_property(problem::IVP{<:AbstractContinuousSystem},
                        property::Property,
                        options::TwoLayerOptions
                       )::Int
    info("Time discretization...")
    Δ = @timing discretize(problem, options[:δ],
        algorithm=options[:discretization], exp_method=options[:exp_method],
        sih_method=options[:sih_method])
    return check_property(Δ, property, options)
end
