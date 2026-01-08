using RecipesBase

#=
The plotting functions in this file are adapted from those in LazySets.jl.
See the relevant docstrings in LazySets for the description of the available options.
=#

using LazySets: plot_recipe,
                PLOT_PRECISION,
                PLOT_POLAR_DIRECTIONS,
                DEFAULT_LABEL,
                DEFAULT_GRID,
                DEFAULT_ASPECT_RATIO,
                DEFAULT_ALPHA,
                DEFAULT_COLOR,
                _extract_limits,
                _extract_extrema,
                _set_auto_limits_to_extrema!,
                _bounding_hyperrectangle,
                _update_plot_limits!,
                isapproxzero,
                _plot_singleton_list

DEFAULT_COLOR_FLOWPIPE = :blue

# heuristics for projecting a reach-set either concretely or lazily to the space
# spanned by vars; the returned set is concrete when the exact projection can be done
# efficiently; in all other cases a lazy set is returned
function _project_reachset(R::AbstractLazyReachSet, vars, ε)
    ST = setrep(R)
    πR = _project_reachset(ST, R, vars)
    if !(πR isa AbstractPolytope) && !isnothing(ε)
        πR = overapproximate(πR, ε)
    end
    return πR
end

# zonotopes are projected concretely unless they have too many generators
function _project_reachset(::Type{<:Zonotope}, R, vars)
    # TODO review/lift the condition on max number of generators after LazySets#2288
    if (ngens(set(R)) <= 15) || (0 ∈ vars)
        # concrete projection is efficient
        πR = set(project(R, vars))

        # zonotopes usually contain lots of redundant generators
        πR = remove_redundant_generators(πR)
    else
        # avoid expensive vertex enumeration
        πR = set(Projection(R, vars))
    end
    return πR
end

function _project_reachset(::Union{Type{<:AbstractZonotope},VPOLY}, R, vars)
    # concrete projection is efficient
    πR = project(R, vars)
    return set(πR)
end

function _project_reachset(::Type{<:AbstractPolyhedron}, R, vars)
    # if the set is polyhedral and two-dimensional
    #   - if time is not required -> don't project
    #   - otherwise, make the projection (lazily)
    if (dim(R) == 2) && vars[1] == 1 && vars[2] == 2
        πR = set(R) # no-op
    else
        πR = set(Projection(R, vars)) # lazy projection
    end
    return πR
end

function _project_reachset(::Type{<:LazySet}, R, vars)
    πR = Projection(R, vars) # lazy projection
    return set(πR)
end

# ---------------------------------------
# Specialization for template reach-sets
# ---------------------------------------
using LazySets: _check_constrained_dimensions

function _is_compatible_projection(R::TemplateReachSet, vars)
    X = set(R)
    for ci in constraints_list(X)
        status = _check_constrained_dimensions(ci, vars)
        status == 0 && return false
    end
    return true
end

# concrete projection for template reach-sets
# TODO specialize for vars = (0, v) and when the direction v is known from a SEV template
function _project_reachset(R::TemplateReachSet, vars)
    # if the exact projection is cheap, we compute it;
    # otherwise we take the lazy projection
    # in both cases 2D bounded sets are eps-close approximated
    if _is_compatible_projection(R, vars)
        πR = project(R, vars)
    else
        πR = Projection(R, vars)
    end
    X = set(πR)
    if dim(X) == 2 && isbounded(X)
        X = overapproximate(X, HPolygon, 1e-3)
    end
    return X
end

# --------------------------------------------
# Specialization for Taylor model reach-sets
# --------------------------------------------

function _project_reachset(T::TaylorModelReachSet, vars, ε)
    R = overapproximate(T, Zonotope)
    return _project_reachset(R, vars, ε)
end

function _check_vars(vars)
    if isnothing(vars)
        throw(ArgumentError("default plotting variables not implemented yet; you need " *
                            "to pass the `vars=(...)` option, e.g. `vars=(0, 1)` to plot variable with " *
                            "index 1 vs. time, or `vars=(1, 2)` to plot variable with index 2 vs. variable with index 1`"))
    end
    D = length(vars)
    @assert (D == 1) || (D == 2) "can only plot in one or two dimensions, " *
                                 "but received $D variable indices where `vars = ` $vars"
end

# ========================
# Reach-set plot recipes
# ========================

@recipe function plot_reachset(R::AbstractReachSet{N};
                               vars=nothing,
                               ε=N(PLOT_PRECISION)) where {N}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR

    return _project_reachset(R, vars, ε)
end

function _plot_reachset_list(list, vars, ε)
    return [_project_reachset(Ri, vars, ε) for Ri in list]
end

# ========================
# Flowpipe plot recipes
# ========================

@recipe function plot_list(list::Union{Flowpipe{N},AbstractVector{<:AbstractReachSet{N}}};
                           vars=nothing, ε=N(PLOT_PRECISION)) where {N}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR_FLOWPIPE
    same_recipe --> true  # LazySets flag; use to show only one label

    @series if !(0 ∈ vars) && (setrep(list) <: AbstractSingleton)
        seriestype --> :scatter
        _plot_singleton_list(list)
    else
        seriestype --> :shape
        _plot_reachset_list(list, vars, ε)
    end
end

# composite flowpipes
@recipe function plot_list(fp::Union{<:HybridFlowpipe{N},<:MixedFlowpipe{N}};
                           vars=nothing, ε=N(PLOT_PRECISION)) where {N}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR_FLOWPIPE
    seriestype --> :shape
    same_recipe --> true  # LazySets flag; use to show only one label

    @series begin
        Xs = LazySet{N}[]
        for F in fp
            append!(Xs, _plot_reachset_list(F, vars, ε))
        end
        Xs
    end
end

# ========================
# Solution plot recipes
# ========================

@recipe function plot_list(sol::ReachSolution{<:Flowpipe{N}};
                           vars=nothing, ε=N(PLOT_PRECISION)) where {N}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR_FLOWPIPE
    same_recipe --> true  # LazySets flag; use to show only one label

    fp = flowpipe(sol)
    @series if !(0 ∈ vars) && (setrep(fp) <: AbstractSingleton)
        seriestype --> :scatter
        _plot_singleton_list(fp)
    else
        seriestype --> :shape
        _plot_reachset_list(fp, vars, ε)
    end
end

# compound solution flowpipes
@recipe function plot_list(sol::Union{ReachSolution{<:MixedFlowpipe{N}},
                                      ReachSolution{<:HybridFlowpipe{N}}};
                           vars=nothing, ε=N(PLOT_PRECISION)) where {N}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR_FLOWPIPE
    seriestype --> :shape
    same_recipe --> true  # LazySets flag; use to show only one label

    @series begin
        Xs = LazySet{N}[]
        for F in flowpipe(sol)
            append!(Xs, _plot_reachset_list(F, vars, ε))
        end
        Xs
    end
end

# TODO new plot recipe to dispatch on ShiftedFlowpipe
# πRi = project(F, i, vars) # project the reach-set

# TODO: refactor with Flowpipe
# TODO extend projection of shifted flowpipes to use lazy projection if needed
@recipe function plot_list(fp::ShiftedFlowpipe{FT,N};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS) where {FT,N}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR_FLOWPIPE
    seriestype --> :shape
    same_recipe --> true  # LazySets flag; use to show only one label

    first = true
    x = Vector{N}()
    y = Vector{N}()
    for i in eachindex(fp.F)
        πRi = project(fp, i, vars) # project the reach-set
        Xi = set(πRi) # extract the set representation

        if Xi isa Intersection
            res = plot_recipe(Xi, ε, Nφ)
        else
            # hard-code overapproximation here to avoid individual
            # compilations for mixed sets
            Pi = isa(Xi, AbstractPolytope) ? Xi : overapproximate(Xi, ε)
            vlist = transpose(hcat(convex_hull(vertices_list(Pi))...))
            if isempty(vlist)
                @warn "overapproximation during plotting was empty"
                continue
            end
            res = vlist[:, 1], vlist[:, 2]
            # add first vertex to "close" the polygon
            push!(res[1], vlist[1, 1])
            push!(res[2], vlist[1, 2])
        end
        if isempty(res)
            continue
        else
            x_new, y_new = res
        end
        if first
            first = false
        else
            push!(x, N(NaN))
            push!(y, N(NaN))
        end
        append!(x, x_new)
        append!(y, y_new)
    end
    return x, y
end

#=
# TODO: default plotting without vars definition
# plot each variable vs time
@recipe function plot_reachset(R::AbstractLazyReachSet{N};
                               ε::N=N(PLOT_PRECISION)
                               ) where {N}
    error("not implemented yet")
end

@recipe function plot_list(list::AbstractVector{<:AbstractReachSet{N}};
                           ε::N=N(PLOT_PRECISION),
                           Nφ::Int=PLOT_POLAR_DIRECTIONS,
                           fast::Bool=true
                          ) where {N}
    error("not implemented yet")
end
=#

# solution from computation without bloating and singleton initial condition
# (eg. with ORBIT) is presened as a scatter plot
@recipe function plot_list(sol::ReachSolution{FT,<:ORBIT};
                           vars=nothing) where {N,RT<:ReachSet{N,<:Singleton{N}},FT<:Flowpipe{N,RT}}
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR_FLOWPIPE
    seriestype --> :scatter
    markershape --> :circle

    _check_vars(vars)
    vx, vy = vars[1], vars[2]
    X(k) = (vx == 0) ? tstart(sol[k]) : element(set(sol[k]))[vx]
    Y(k) = (vy == 0) ? tstart(sol[k]) : element(set(sol[k]))[vy]

    x = Vector{N}()
    y = Vector{N}()
    for k in eachindex(sol)
        x_new = X(k)
        y_new = Y(k)
        append!(x, x_new)
        append!(y, y_new)
    end
    return x, y
end
