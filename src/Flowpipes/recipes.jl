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

import LazySets: _plot_singleton_list_1D,
                 _plot_singleton_list_2D

# heuristics for projecting a reach-set either concretely or lazily to the space
# spanned by vars; the returned set is concrete when the exact projection can be done
# efficiently; in all other cases a lazy set is returned
function _project_reachset(R::AbstractLazyReachSet, vars)
    ST = setrep(R)

    if (ST <: AbstractZonotope) || (ST <: VPolytope) || (ST <: VPolygon)
        # concrete projection is efficient
        πR = project(R, vars)
        X = set(πR)

    elseif (ST <: AbstractPolyhedron) && (dim(R) == 2)
        # if the set is polyhedral and two-dimensional
        #   - if time is not requried -> don't project
        #   - otherwise, make the projection (lazily)
        X = 0 ∈ vars ? set(Projection(R, vars)) : set(R)

    else
        πR = Projection(R, vars) # lazy projection
        X = set(πR)
    end
    return X
end

function _project_reachset(T::TaylorModelReachSet, vars)
    R = overapproximate(T, Zonotope)
    _project_reachset(R, vars)
end

# templates have their own projection heuristics
function _project_reachset(R::TemplateReachSet, vars)
    πR = project(R, vars)
    X = set(πR)
end

function _check_vars(vars)
    if vars == nothing
        throw(ArgumentError("default ploting variables not implemented yet; you need " *
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

@recipe function plot_reachset(R::AbstractLazyReachSet{N};
                               vars=nothing,
                               ε=N(PLOT_PRECISION)) where {N}
    _check_vars(vars)
    X = _project_reachset(R, vars)

    dim(X) == 1 && return plot_recipe(X, ε)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR

    # extract limits and extrema of already plotted sets
    p = plotattributes[:plot_object]
    lims = _extract_limits(p)
    extr = _extract_extrema(p)

    if !isbounded(X)
        _set_auto_limits_to_extrema!(lims, extr)
        X = intersection(X, _bounding_hyperrectangle(lims, eltype(X)))

    elseif length(p) > 0
        # if there is already a plotted set and the limits are fixed,
        # automatically adjust the axis limits (e.g. after plotting a unbounded set)
        _update_plot_limits!(lims, X)
    end

    xlims --> lims[:x]
    ylims --> lims[:y]

    x, y = plot_recipe(X, ε)
    if !isempty(x)
        m = length(x)
        if m == 1
            seriestype := :scatter
        elseif m == 2 && isapproxzero((x[1] - x[2])^2 + (y[1] - y[2])^2)
            seriestype := :scatter
        else
            seriestype := :shape
        end
    end
    return (x, y)
end

# TODO use LazySets._plot_singleton_list_1D if the reach-sets are stored as a struct
# array. Otherwise, we could use pass the list through x -> set(x)
function _plot_singleton_list_1D(list::Union{Flowpipe{N}, AbstractVector{RN}}) where {N, RN<:AbstractReachSet{N}}
    m = length(list)

    x = Vector{N}(undef, m)
    y = zeros(N, m)

    @inbounds for (i, Xi) in enumerate(list)
        p = element(set(Xi))
        x[i] = p[1]
    end
    x, y
end

function _plot_singleton_list_2D(list::Union{Flowpipe{N}, AbstractVector{RN}}) where {N, RN<:AbstractReachSet{N}}
    m = length(list)
    x = Vector{N}(undef, m)
    y = Vector{N}(undef, m)

    @inbounds for (i, Xi) in enumerate(list)
        p = element(set(Xi))
        x[i] = p[1]
        y[i] = p[2]
    end
    x, y
end

function _plot_reachset_list(list, N, vars, ε, Nφ)
    first = true
    x = Vector{N}()
    y = Vector{N}()
    for Ri in list
        Xi = _project_reachset(Ri, vars)

        if Xi isa Intersection
            xcoords, ycoords = plot_recipe(Xi, ε, Nφ)

        else
            Pi = isconcretetype(typeof(Xi)) ? Xi : overapproximate(Xi, ε)
            vlist = convex_hull(vertices_list(Pi))
            m = length(vlist)
            if m == 0
                @warn "overapproximation during plotting was empty"
                continue
            end
            xcoords = Vector{N}(undef, m)
            ycoords = Vector{N}(undef, m)
            @inbounds for (i, v) in enumerate(vlist)
                xcoords[i] = v[1]
                ycoords[i] = v[2]
            end
            if m > 1
                # add first vertex to "close" the polygon
                push!(xcoords, xcoords[1])
                push!(ycoords, ycoords[1])
            end
        end
        isempty(xcoords) && continue

        x_new = xcoords
        y_new = ycoords

        if first
            first = false
        else
            push!(x, N(NaN))
            push!(y, N(NaN))
        end
        append!(x, x_new)
        append!(y, y_new)
    end
    x, y
end

# ========================
# Flowpipe plot recipes
# ========================

@recipe function plot_list(list::Union{Flowpipe{N}, AbstractVector{RN}};
                           vars=nothing,
                           ε=N(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS) where {N, RN<:AbstractReachSet{N}}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR

    if !(0 ∈ vars) && (setrep(list) <: AbstractSingleton)
        seriestype --> :scatter
        _plot_singleton_list(list)
    else
        seriestype --> :shape
        _plot_reachset_list(list, N, vars, ε, Nφ)
    end
end

# compound flowpipes
@recipe function plot_list(fp::Union{HF, MF};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS
                          ) where {N, HF<:HybridFlowpipe{N}, MF<:MixedFlowpipe{N}}

    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype --> :shape

    x = Vector{N}()
    y = Vector{N}()

    for F in fp
        x_new, y_new = _plot_reachset_list(F, N, vars, ε, Nφ)
        append!(x, x_new)
        append!(y, y_new)
    end
    x, y
end


# ========================
# Solution plot recipes
# ========================

@recipe function plot_list(sol::ReachSolution{FT};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS
                          ) where {N, FT<:Flowpipe{N}}

    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR

    fp = flowpipe(sol)
    if !(0 ∈ vars) && (setrep(fp) <: AbstractSingleton)
        seriestype --> :scatter
        _plot_singleton_list(fp)
    else
        seriestype --> :shape
        _plot_reachset_list(fp, N, vars, ε, Nφ)
    end
end

# compound solution flowpipes
@recipe function plot_list(sol::Union{SMF, SHF};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS
                          ) where {N, MF<:MixedFlowpipe{N}, HF<:HybridFlowpipe{N}, SMF<:ReachSolution{MF}, SHF<:ReachSolution{HF}}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype --> :shape

    first = true
    x = Vector{N}()
    y = Vector{N}()
    for F in flowpipe(sol)
    for (i, Ri) in enumerate(F)
        if isa(F, ShiftedFlowpipe)
            # TODO refactor; this is needed to support ShiftedFlowpipe
            πRi = project(F, i, vars) # project the reach-set
            Xi = set(πRi) # extract the set representation
        else
            Xi = _project_reachset(Ri, vars)
        end

        if Xi isa Intersection
            res = plot_recipe(Xi, ε, Nφ)
        else
            # hard-code overapproximation here to avoid individual
            # compilations for mixed sets
            Pi = isa(Xi, AbstractPolygon) ? Xi : overapproximate(Xi, ε)
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
    end
    x, y
end

# TODO: refactor with Flowpipe
# TODO extend projection of shifted flowpipes to use lazy projection if needed
@recipe function plot_list(fp::ShiftedFlowpipe{N};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS
                          ) where {N}
    _check_vars(vars)

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype --> :shape

    first = true
    x = Vector{N}()
    y = Vector{N}()
    for (i, Ri) in enumerate(fp)
        πRi = project(fp, i, vars) # project the reach-set
        Xi = set(πRi) # extract the set representation

        if Xi isa Intersection
            res = plot_recipe(Xi, ε, Nφ)
        else
            # hard-code overapproximation here to avoid individual
            # compilations for mixed sets
            Pi = isa(Xi, AbstractPolygon) ? Xi : overapproximate(Xi, ε)
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
    x, y
end

#=
# TODO: default plotting without vars definition
# plot each variable vs time
@recipe function plot_reachset(R::AbstractLazyReachSet{N};
                               ε::N=N(PLOT_PRECISION)
                               ) where {N<:Real, D, M<:Integer}
    error("not implemented yet")
end

@recipe function plot_list(list::AbstractVector{RN};
                           ε::N=N(PLOT_PRECISION),
                           Nφ::Int=PLOT_POLAR_DIRECTIONS,
                           fast::Bool=true
                          ) where {N<:Real, RN<:AbstractReachSet{N}}
    error("not implemented yet")
end
=#

# solution from computation without bloating and singleton initial condition
# (eg. with ORBIT) is presened as a scatter plot
@recipe function plot_list(sol::ReachSolution{FT, <:ORBIT};
                           vars=nothing) where {N, RT<:ReachSet{N, <:Singleton{N}}, FT<:Flowpipe{N, RT}}
   label --> DEFAULT_LABEL
   grid --> DEFAULT_GRID
   if DEFAULT_ASPECT_RATIO != :none
       aspect_ratio --> DEFAULT_ASPECT_RATIO
   end
   seriesalpha --> DEFAULT_ALPHA
   seriescolor --> DEFAULT_COLOR
   seriestype --> :scatter
   markershape --> :circle

   _check_vars(vars)
   vx, vy = vars[1], vars[2]
   X(k) = (vx == 0) ? tstart(sol[k]) : element(set(sol[k]))[vx]
   Y(k) = (vy == 0) ? tstart(sol[k]) : element(set(sol[k]))[vy]

   x = Vector{N}()
   y = Vector{N}()
   for k in 1:length(sol)
       x_new = X(k)
       y_new = Y(k)
       append!(x, x_new)
       append!(y, y_new)
   end
   x, y
end
