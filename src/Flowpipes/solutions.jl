# ================================
# Abstract types
# ================================

"""
    AbstractSolution

Abstract supertype of all solution types of a rechability problem.
"""
abstract type AbstractSolution end

# ================================
# Property checking problem
# ================================

"""
    CheckSolution{T} <: AbstractSolution

Type that wraps the solution of a verification problem.

### Fields

TODO

### Notes

This type contains the answer if the property is satisfied, and if not, it
contains the index at which the property might be violated for the first time.
"""
struct CheckSolution{PT, ST<:AbstractPost} <: AbstractSolution
    property::PT
    satisfied::Bool
    vidx::Int
    vtspan::TimeInterval
    solver::ST
end

property(sol::CheckSolution) = sol.property
satisfied(sol::CheckSolution) = sol.satisfied
violation_index(sol::CheckSolution) = sol.vidx
violation_tspan(sol::CheckSolution) = sol.vtspan
solver(sol::CheckSolution) = sol.solver

# ================================
# Reachability problem
# ================================

"""
    ReachSolution{FT<:AbstractFlowpipe, ST<:AbstractPost} <: AbstractSolution

Type that wraps the solution of a reachability problem as a sequence of lazy
sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of [`AbstractReachSet`](@ref)s
- `options`  -- the dictionary of options
"""
struct ReachSolution{FT<:AbstractFlowpipe, ST<:AbstractPost} <: AbstractSolution
    F::FT
    solver::ST
    ext::Dict{Symbol, Any} # dictionary used by extensions
end

# constructor from empty extension dictionary
function ReachSolution(F::FT, solver::ST) where {FT<:AbstractFlowpipe, ST<:AbstractPost}
    return ReachSolution(F, solver, Dict{Symbol, Any}())
end

# getter functions
flowpipe(sol::ReachSolution) = sol.F
tstart(sol::ReachSolution) = tstart(sol.F)
tend(sol::ReachSolution) = tend(sol.F)
tspan(sol::ReachSolution) = tspan(sol.F)
LazySets.dim(sol::ReachSolution) = dim(sol.F) # TODO: keep for hybrid?

# iteration and indexing iterator interface
array(sol::ReachSolution) = array(sol.F)
Base.iterate(sol::ReachSolution) = iterate(sol.F)
Base.iterate(sol::ReachSolution, state) = iterate(sol.F, state)
Base.length(sol::ReachSolution) = length(sol.F)
Base.first(sol::ReachSolution) = first(sol.F)
Base.last(sol::ReachSolution) = last(sol.F)
Base.firstindex(sol::ReachSolution) = 1
Base.lastindex(sol::ReachSolution) = length(sol.F)
Base.getindex(sol::ReachSolution, i::Int) = getindex(sol.F, i)
Base.getindex(sol::ReachSolution, i::Number) = getindex(sol.F, i)
Base.getindex(sol::ReachSolution, I::AbstractVector) = getindex(sol.F, I)

# evaluation interface
Base.getindex(sol::ReachSolution, t::Float64) = getindex(sol.F, t)
(sol::ReachSolution)(t::Float64) = sol.F(t)
(sol::ReachSolution)(t::Number) = sol.F(t)
(sol::ReachSolution)(dt::IA.Interval{Float64}) = sol.F(dt)

#=
function project(rs::ReachSolution, M::AbstractMatrix)
    Yk = [project(X, M) for X in rs.Xk]
    return ReachSolution(Yk, rs.options)
end
=#

#=

"""
    plot_sol(sol::ReachSolution; ...)

Plots the solution of a reachability problem in 2D with the given options.

### Input

- `sol`            --  the solution of a reachability problem, projected into
                       two dimensions
- `seriescolor`     -- (optional, default: "blue"): color of the polygons; by
                       default, the same color is used for all of them
- `label`           -- (optional, default: nothing): the legend label
- `grid`            -- (optional, default: true): to use gridlines or not
- `alpha`           -- (optional, default: 0.5): the transparency of the polygons
- `indices`         -- (optional, default: nothing): if given, plot only a sample of
                       the array; this option can be passed through `sol.options`
- `vars`            -- (optional, default: nothing): if given, label the x and y
                       axes with the variables being plotted; this option can be
                       passed through `sol.options`, or as a pair of integers,
                       where 0 indicates the time variable
- `use_subindices`  -- (optional, default: `false`) if `false`, use subindices
                       for the labels, e.g. `x1` is displayed as `x₁`

### Notes

To define your own x and y labels, use the `xlabel` (resp. `ylabel`) keyword
argument. For additional options, consult the Plots.jl reference manual.

# TODO: add warning if too many sets are tried to be plotted, and use sampling
# if it is the case.
"""
@recipe function plot_sol(sol::ReachSolution;
                          seriescolor=:auto,
                          fillcolor=:auto,
                          seriestype=:shape,
                          label="", grid=true, alpha=0.5,
                          indices=nothing, vars=nothing,
                          use_subindices=false)
    @assert dim(set(sol.Xk[1])) == 2 "we only support plotting 2D sets"

    options = check_aliases_and_add_default_value(sol)

    if vars != nothing
        vars = add_plot_labels(vars, use_subindices=use_subindices)
        xguide --> vars[1]; yguide --> vars[2]
    elseif options[:plot_vars] != nothing
        vars = add_plot_labels(options[:plot_vars], use_subindices=use_subindices)
        xguide --> vars[1]; yguide --> vars[2]
    end

    if indices == nothing
        indices = options[:plot_indices]
    end

    N = length(indices)
    # rule of thumb for linecolor and fillcolor overlap when plotting many sets
    if N < 300
        linecolor --> :black
    else
        linecolor --> :match
    end

    # Using single list and NaN separators
    vlist = Vector{Vector{Float64}}()
    for i in indices
        append!(vlist, convex_hull(vertices_list(set(sol.Xk[i]))))
        push!(vlist, [NaN; NaN])
    end
    vlist = hcat(vlist...)'
    vlist[:, 1], vlist[:, 2]

end


"""
    ReachSolution{SN, RSN<:AbstractReachSet{SN}} <: AbstractSolution

Type that wraps the solution of a reachability problem as a sequence of lazy
sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of [`AbstractReachSet`](@ref)s
- `options`  -- the dictionary of options
"""
struct ReachSolution{SN, RSN<:AbstractReachSet{SN}, T} <: AbstractSolution
    Xk::Flowpipe
    solver::T
end

# constructor with no options
ReachSolution(Xk::Vector{RSN}) where {SN, RSN<:AbstractReachSet{SN}} =
    ReachSolution{SN, RSN}(Xk, Options())

function project(rs::ReachSolution, M::AbstractMatrix)
    Yk = [project(X, M) for X in rs.Xk]
    return ReachSolution(Yk, rs.options)
end

"""
    plot_sol(sol::ReachSolution; ...)

Plots the solution of a reachability problem in 2D with the given options.

### Input

- `sol`            --  the solution of a reachability problem, projected into
                       two dimensions
- `seriescolor`     -- (optional, default: "blue"): color of the polygons; by
                       default, the same color is used for all of them
- `label`           -- (optional, default: nothing): the legend label
- `grid`            -- (optional, default: true): to use gridlines or not
- `alpha`           -- (optional, default: 0.5): the transparency of the polygons
- `indices`         -- (optional, default: nothing): if given, plot only a sample of
                       the array; this option can be passed through `sol.options`
- `vars`            -- (optional, default: nothing): if given, label the x and y
                       axes with the variables being plotted; this option can be
                       passed through `sol.options`, or as a pair of integers,
                       where 0 indicates the time variable
- `use_subindices`  -- (optional, default: `false`) if `false`, use subindices
                       for the labels, e.g. `x1` is displayed as `x₁`

### Notes

To define your own x and y labels, use the `xlabel` (resp. `ylabel`) keyword
argument. For additional options, consult the Plots.jl reference manual.
"""
@recipe function plot_sol(sol::ReachSolution;
                          seriescolor=:auto,
                          fillcolor=:auto,
                          seriestype=:shape,
                          label="", grid=true, alpha=0.5,
                          indices=nothing, vars=nothing,
                          use_subindices=false)
    @assert dim(set(sol.Xk[1])) == 2 "we only support plotting 2D sets"

    options = check_aliases_and_add_default_value(sol)

    if vars != nothing
        vars = add_plot_labels(vars, use_subindices=use_subindices)
        xguide --> vars[1]; yguide --> vars[2]
    elseif options[:plot_vars] != nothing
        vars = add_plot_labels(options[:plot_vars], use_subindices=use_subindices)
        xguide --> vars[1]; yguide --> vars[2]
    end

    if indices == nothing
        indices = options[:plot_indices]
    end

    N = length(indices)
    # rule of thumb for linecolor and fillcolor overlap when plotting many sets
    if N < 300
        linecolor --> :black
    else
        linecolor --> :match
    end

    # Using single list and NaN separators
    vlist = Vector{Vector{Float64}}()
    for i in indices
        append!(vlist, convex_hull(vertices_list(set(sol.Xk[i]))))
        push!(vlist, [NaN; NaN])
    end
    vlist = hcat(vlist...)'
    vlist[:, 1], vlist[:, 2]

end
=#
