using ..Utils: LDS, CLCDS

using LazySets: CachedMinkowskiSumArray,
                 isdisjoint

import LazySets.Approximations: overapproximate

"""
combine_cpas(cpa1::CartesianProductArray{N, S}, cpa2::CartesianProductArray{N, S},
                       blocks1::Vector{Int}, blocks2::Vector{Int}) where {N, S<:LazySet{N}}

Compose high-dimensional cartesian product array from 2 low-dimensional CPA's

### Input

- `cpa1`   -- Cartesian Product Array
- `cpa2`   -- Cartesian Product Array
- `blocks1`   -- Block structure of cpa1 in sense of high-dimensional set
- `blocks2`   -- Block structure of cpa2 in sense of high-dimensional set

### Output

High-dimensional cartesian product array.
"""
function combine_cpas(cpa1::CartesianProductArray{N, S}, cpa2::CartesianProductArray{N, S},
                       blocks1::Vector{Int}, blocks2::Vector{Int}) where {N, S<:LazySet{N}}
    result = Vector{S}(undef, length(array(cpa1)) + length(array(cpa2)))
    result[blocks1] = array(cpa1)
    result[blocks2] = array(cpa2)
    return CartesianProductArray(result)
end

"""
    reach(problem, options)

Interface to reachability algorithms for an LTI system.

### Input

- `problem`   -- initial value problem
- `options`   -- additional options

### Output

A sequence of [`SparseReachSet`](@ref)s.

### Notes

The numeric type of the system's coefficients and the set of initial states
is inferred from the first parameter of the system (resp. lazy set), ie.
`NUM = first(typeof(problem.s).parameters)`.
"""
function reach_mixed(problem::Union{IVP{<:CLDS{NUM}, <:LazySet{NUM}},
                              IVP{<:CLCDS{NUM}, <:LazySet{NUM}}},
               options::TwoLayerOptions
              ) where {NUM <: Real}
    # optional matrix conversion
    problem = matrix_conversion(problem, options)

    # list containing the arguments passed to any reachability function
    args = []

    # coefficients matrix
    A = problem.s.A
    push!(args, A)

    # determine analysis mode (sparse/dense) for lazy_expm mode
    if A isa SparseMatrixExp
        push!(args, Val(options[:assume_sparse]))
    end

    n = statedim(problem)
    blocks = options[:blocks]
    partition = convert_partition(options[:partition])
    T = options[:T]
    N = (T == Inf) ? nothing : ceil(Int, T / options[:δ])
    # Cartesian decomposition of the initial set
    if length(partition) == 1 && length(partition[1]) == n &&
            options[:block_options_init] == LinearMap
        info("- Skipping decomposition of X0")
        Xhat0 = LazySet{NUM}[problem.x0]
    else
        info("- Decomposing X0")
        @timing begin
            Xhat0 = array(decompose(problem.x0, options[:partition],
                                    options[:block_options_init]))
        end
    end

    # compute dimensions
    dimensions = compute_dimensions(partition, blocks)

    # determine output function: linear map with the given matrix
    output_function = options[:output_function] != nothing ?
        (x -> options[:output_function] * x) :
        nothing

    # preallocate output vector
    if output_function == nothing
        res_type = SparseReachSet{CartesianProductArray{NUM, LazySet{NUM}}}
    else
        res_type = SparseReachSet{Hyperrectangle{NUM}}
    end
    res = (N == nothing) ? Vector{res_type}() : Vector{res_type}(undef, N)

    # shortcut if only the initial set is required
    if N == 1
        res[1] = res_type(
            CartesianProductArray{NUM, LazySet{NUM}}(Xhat0[blocks]),
            zero(NUM), options[:δ])
        return res
    end
    push!(args, Xhat0)

    # inputs
    if !options[:assume_homogeneous] && inputdim(problem) > 0
        U = inputset(problem)
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
    push!(args, overapproximate_fun)

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

    # output function
    push!(args, output_function)

    # add mode-specific block(s) argument
    push!(args, blocks)
    push!(args, partition)
    push!(args, dimensions)

    # time step
    push!(args, options[:δ])

    # termination function
    invariant = stateset(problem.s)
    vars = [e for block in blocks for e in partition[block]]
    if invariant isa Universe
        invariant = Universe(length(vars))
    else
        invariant = LazySets.Approximations.project(invariant, vars)
    end
    termination = get_termination_function_out(N, invariant)
    if length(options[:guards_proj]) > 0
        push!(args, options[:guards_proj])
    else
        push!(args, [Universe{NUM}(length(vars))])
    end
    push!(args, options[:block_options_iter])
    push!(args, vars)
    push!(args, termination)

    info("- Computing successors")

    # progress meter
    progress_meter = (N != nothing) ?
        Progress(N, 1, "Computing successors ") :
        nothing
    push!(args, progress_meter)

    # preallocated result
    push!(args, res)

    # call the adequate function with the given arguments list
    @timing begin
        index, skip = reach_blocks!(args...)
    end

    # shrink result array
    if skip || index < N
        if N != nothing
            info("termination after only $index of $N steps")
        end
        deleteat!(res, (skip ? index : index + 1):length(res))
    end

    # return the result
    return res
end

function reach_mixed(problem::Union{IVP{<:CLCS{NUM}, <:LazySet{NUM}},
                              IVP{<:CLCCS{NUM}, <:LazySet{NUM}}},
               options::TwoLayerOptions
              ) where {NUM <: Real}
    info("Time discretization...")
    Δ = @timing discretize(problem, options[:δ],
        algorithm=options[:discretization], exp_method=options[:exp_method],
        sih_method=options[:sih_method])
    return reach_mixed(Δ, options)
end

function get_termination_function_out(N::Nothing, invariant::Universe)
    termination_function = (k, set, t0) -> termination_unconditioned(set)
    warn("no termination condition specified; the reachability analysis will " *
         "not terminate")
    return termination_function
end

function get_termination_function_out(N::Int, invariant::Universe)
    return (k, set, t0) -> termination_N_out(N, k, t0, set)
end

function get_termination_function_out(N::Nothing, invariant::LazySet)
    return (k, set, t0) -> termination_inv_out(invariant, set, t0)
end

function get_termination_function_out(N::Int, invariant::LazySet)
    return (k, set, t0) -> termination_inv_N_out(N, invariant, k, set, t0)
end

function termination_unconditioned_out(set)
    return (false, false, set)
end

function termination_N_out(N, k, t0, set)
    return (k >= N, false, set)
end

function termination_inv_out(inv, set, t0)
    l_inters = Intersection(set, inv)
    if isempty(l_inters)
        return (true, true, EmptySet())
    else
        return (false, false, l_inters)
    end
end

function termination_inv_N_out(N, inv, k, set, t0)
    if k >= N
        return (true, false, EmptySet())
    end
    l_inters = Intersection(set, inv)
    if isempty(l_inters)
        return (true, true, EmptySet())
    else
        return (false, false, l_inters)
    end
end
