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
                _update_plot_limits!

# heuristics for projecting a reach-set either concretely or lazily
function _project_reachset(R::AbstractLazyReachSet{N}, vars, ε=N(PLOT_PRECISION)) where {N}
    if setrep(R) <: AbstractZonotope
        # concrete projection is efficient for zonotopic sets
        πR = project(R, vars)
        X = set(πR)
    elseif (setrep(R) <: AbstractPolyhedron) && (dim(R) == 2)

        # TODO : << cleanup plotting heuristics for polyhedral sets with / without projection

        # 2D polyhedral sets do not need to be projected unless one of the
        # coordinates of interest is time; in that case, we take the lazy projection
        # and then overapproximate with a box without loss
        # NOTE: ε option is currently ignored
        #if vars == 1:2
        #    X = set(R)
        #else
        #    πR = Projection(R, vars)
        #    X = overapproximate(set(πR), Hyperrectangle)
        #end

        # - if the set is polyhedral and 2D
        #   - if time is not requried -> don't project
        #   - otherwise, make the projection (lazily)
        if 0 ∈ vars
            πR = Projection(R, vars) # lazy projection
            #X = overapproximate(set(πR), ε)
            X = set(πR)
        else
            X = set(R)
        end
    else
        πR = Projection(R, vars) # lazy projection
        X = overapproximate(set(πR), ε)
    end
    return X
end

function _project_reachset(T::TaylorModelReachSet, vars, ε=N(PLOT_PRECISION))
    R = overapproximate(T, Zonotope)
    _project_reachset(R, vars, ε)
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

# The type annotation NTuple in vars is removed because of this warning:
# ┌ Warning: Type annotations on keyword arguments not currently supported in recipes. Type information has been discarded
# └ @ RecipesBase ~/.julia/packages/RecipesBase/zBoFG/src/RecipesBase.jl:112
@recipe function plot_reachset(R::AbstractLazyReachSet{N};
                               vars=nothing,
                               ε=N(PLOT_PRECISION)
                               ) where {N<:Real}


    _check_vars(vars)
    X = _project_reachset(R, vars, ε)

    if dim(X) == 1
        plot_recipe(X, ε)
    else
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

        # if there is already a plotted set and the limits are fixed,
        # automatically adjust the axis limits (e.g. after plotting a unbounded set)
        elseif length(p) > 0
            _update_plot_limits!(lims, X)
        end

        xlims --> lims[:x]
        ylims --> lims[:y]

        res = plot_recipe(X, ε)
        if isempty(res)
            res
        else
            x, y = res
            if length(x) == 1 || norm([x[1], y[1]] - [x[2], y[2]]) ≈ 0
                seriestype := :scatter
            else
                seriestype := :shape
            end
            x, y
        end
    end
end

@recipe function plot_list(list::AbstractVector{RN};
                           vars=nothing,
                           ε=N(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS,
                           fast=true
                          ) where {N<:Real, RN<:AbstractReachSet{N}}

    _check_vars(vars)

    if fast
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
        for Ri in list
            Xi = _project_reachset(Ri, vars, ε)

            if Xi isa Intersection
                res = plot_recipe(Xi, ε, Nφ)
            else
                # hard-code overapproximation here to avoid individual
                # compilations for mixed sets
                Pi = overapproximate(Xi, ε)
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
    else
        for Ri in list
            Xi = _project_reachset(Ri, vars, ε)
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
    end
end

# This function is from LazySets.jl. See the docstring in LazySets for the description
# of the available options.
@recipe function plot_list(fp::Flowpipe{N};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS,
                           fast=true
                          ) where {N}

    _check_vars(vars)

    if fast
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
        for Ri in fp
            Xi = _project_reachset(Ri, vars, ε)

            if Xi isa Intersection
                res = plot_recipe(Xi, ε, Nφ)
            else
                # hard-code overapproximation here to avoid individual
                # compilations for mixed sets
                Pi = overapproximate(Xi, ε)
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
    else
        for Ri in fp
            Xi = _project_reachset(Ri, vars, ε)
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
    end
end

# compound flowpipes
@recipe function plot_list(fp::Union{HF, MF};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS,
                           fast=true
                          ) where {N, HF<:HybridFlowpipe{N}, MF<:MixedFlowpipe{N}}

    _check_vars(vars)

    if fast
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
        for F in fp
        for Ri in F
            Xi = _project_reachset(Ri, vars, ε)

            if Xi isa Intersection
                res = plot_recipe(Xi, ε, Nφ)
            else
                # hard-code overapproximation here to avoid individual
                # compilations for mixed sets
                Pi = overapproximate(Xi, ε)
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
    else
        for F in fp
        for Ri in F
            Xi = _project_reachset(Ri, vars, ε)
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
        end
    end
end

@recipe function plot_list(sol::ReachSolution{FT};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS,
                           fast=true
                          ) where {N, FT<:Flowpipe{N}}

    _check_vars(vars)
    if fast
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
        for Ri in flowpipe(sol)
            Xi = _project_reachset(Ri, vars, ε)
            if Xi isa Intersection
                res = plot_recipe(Xi, ε, Nφ)
            else
                # hard-code overapproximation here to avoid individual
                # compilations for mixed sets
                Pi = overapproximate(Xi, ε)
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
    else
        for Ri in flowpipe(sol)
            Xi = _project_reachset(Ri, vars, ε)
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
    end
end

# compound solution flowpipes
@recipe function plot_list(sol::Union{SMF, SHF};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS,
                           fast=true
                          ) where {N, MF<:MixedFlowpipe{N}, HF<:HybridFlowpipe{N}, SMF<:ReachSolution{MF}, SHF<:ReachSolution{HF}}
    _check_vars(vars)
    if fast
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
                Xi = _project_reachset(Ri, vars, ε)
            end

            if Xi isa Intersection
                res = plot_recipe(Xi, ε, Nφ)
            else
                # hard-code overapproximation here to avoid individual
                # compilations for mixed sets
                Pi = overapproximate(Xi, ε)
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
    else
        for F in flowpipe(sol)
        for (i, Ri) in enumerate(F)
            if isa(F, ShiftedFlowpipe)
                # TODO refactor; this is needed to support ShiftedFlowpipe
                πRi = project(F, i, vars) # project the reach-set
                Xi = set(πRi) # extract the set representation
            else
                Xi = _project_reachset(Ri, vars, ε)
            end
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
        end
    end
end

# TODO: refactor with Flowpipe
# TODO extend projection of shifted flowpipes to use lazy projection if needed
@recipe function plot_list(fp::ShiftedFlowpipe{N};
                           vars=nothing,
                           ε=Float64(PLOT_PRECISION),
                           Nφ=PLOT_POLAR_DIRECTIONS,
                           fast=true
                          ) where {N}
    _check_vars(vars)
    if fast
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
                Pi = overapproximate(Xi, ε)
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
    else
        for (i, Ri) in enumerate(fp)
            πRi = project(fp, i, vars) # project the reach-set
            Xi = set(πRi) # extract the set representation
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
    end
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
@recipe function plot_list(sol::ReachSolution{FT, <:ORBIT}; vars=nothing) where {N, RT<:ReachSet{N, <:Singleton{N}}, FT<:Flowpipe{N, RT}}
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
