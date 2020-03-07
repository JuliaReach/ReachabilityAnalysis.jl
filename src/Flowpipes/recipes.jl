using RecipesBase

using LazySets: plot_recipe,
                PLOT_PRECISION,
                DEFAULT_LABEL,
                DEFAULT_GRID,
                DEFAULT_ASPECT_RATIO,
                DEFAULT_ALPHA,
                DEFAULT_COLOR,
                _extract_limits,
                _extract_extrema,
                _set_auto_limits_to_extrema!,
                _bounding_hyperrectangle

# This function is from LazySets.jl. See the docstring in LazySets for the description
# of the available optoions.
@recipe function plot_reachset(R::AbstractLazyReachSet{N}, ε::N=N(PLOT_PRECISION)) where {N<:Real}
    X = set(R)
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

# This function is from LazySets.jl. See the docstring in LazySets for the description
# of the available optoions.
@recipe function plot_list(list::AbstractVector{RN}, ε::N=N(PLOT_PRECISION),
                           Nφ::Int=PLOT_POLAR_DIRECTIONS, fast::Bool=false
                          ) where {N<:Real, RN<:AbstractReachSet{N}}
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
            Xi = set(Ri)
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
            Xi = set(Ri)
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
    end
end
