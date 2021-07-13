# ================================
# Flowpipe
# ================================

"""
    Flowpipe{N, RT<:AbstractReachSet{N}, VRT<:AbstractVector{RT}} <: AbstractFlowpipe

Type that wraps a flowpipe, which is an iterable collection of reach-sets that
behaves like their set union.

### Fields

- `Xk`  -- array of reach-sets
- `ext` -- extension dictionary; field used by extensions

### Notes

The dimension of the flowpipe corresponds to the dimension of the underlying
reach-sets; in this type, it is is assumed that the dimension is the same for
the different reach-sets.
"""
struct Flowpipe{N, RT<:AbstractReachSet{N}, VRT<:AbstractVector{RT}} <: AbstractFlowpipe
    Xk::VRT
    ext::Dict{Symbol, Any}
end

# getter functions
@inline array(fp::Flowpipe) = fp.Xk
@inline flowpipe(fp::Flowpipe) = fp

# constructor from empty extension dictionary
function Flowpipe(Xk::AbstractVector{RT}) where {N, RT<:AbstractReachSet{N}}
    return Flowpipe(Xk, Dict{Symbol, Any}())
end

# undef initializer given a set type
function Flowpipe(::UndefInitializer, ST::Type{<:LazySet{N}}, k::Int) where {N}
    return Flowpipe(Vector{ReachSet{N, ST}}(undef, k), Dict{Symbol, Any}())
end

# undef initializer given a reach-set type
function Flowpipe(::UndefInitializer, RT::Type{<:AbstractReachSet}, k::Int)
    return Flowpipe(Vector{RT}(undef, k), Dict{Symbol, Any}())
end

# undef initializer given a continuous post-operator
function Flowpipe(::UndefInitializer, cpost::AbstractContinuousPost, k::Int)
    RT = rsetrep(cpost)
    return Flowpipe(Vector{RT}(undef, k), Dict{Symbol, Any}())
end

# constructor from a single reach-set
Flowpipe(R::AbstractReachSet) = Flowpipe([R])

function Base.similar(fp::Flowpipe{N, RT, VRT}) where {N, RT, VRT}
   return Flowpipe(VRT())
end

Base.IndexStyle(::Type{<:Flowpipe}) = IndexLinear()
Base.eltype(::Flowpipe{N, RT}) where {N, RT} = RT
Base.size(fp::Flowpipe) = (length(fp.Xk),)
Base.view(fp::Flowpipe, args...) = view(fp.Xk, args...)
Base.push!(fp::Flowpipe, args...) = push!(fp.Xk, args...)
Base.keys(fp::Flowpipe) = keys(fp.Xk)

numtype(::Flowpipe{N}) where {N} = N
setrep(fp::Flowpipe{N, RT}) where {N, RT} = setrep(RT)
setrep(::Type{<:Flowpipe{N, RT}}) where {N, RT} = setrep(RT)
rsetrep(fp::Flowpipe{N, RT}) where {N, RT} = RT
rsetrep(::Type{<:Flowpipe{N, RT}}) where {N, RT} = RT
numrsets(fp::Flowpipe) = length(fp)

# getter functions for hybrid systems

"""
    location(F::Flowpipe)

Return the location of a flowpipe within a hybrid system, or `missing` if it is
not defined.

### Input

- `F` -- flowpipe

### Output

The `:loc_id` value of stored in the flowpipe's extension field.
"""
function location(fp::Flowpipe)
    return get(fp.ext, :loc_id, missing)
end

# evaluate a flowpipe at a given time point: gives a reach set
# here it would be useful to layout the times contiguously in a vector
# (see again array of struct vs struct of array)
function (fp::Flowpipe)(t::Number)
    Xk = array(fp)
    @inbounds for (i, X) in enumerate(Xk)
        if t ∈ tspan(X) # exit on the first occurrence
            if i < length(Xk) && t ∈ tspan(Xk[i+1])
                return view(Xk, i:i+1)
            else
                return X
            end
        end
    end
    throw(ArgumentError("time $t does not belong to the time span, " *
                        "$(tspan(fp)), of the given flowpipe"))
end

# evaluate a flowpipe at a given time interval: gives possibly more than one reach set
# i.e. first and last sets and those in between them
function (fp::Flowpipe)(dt::TimeInterval)
    idx = Vector{Int}()
    α = inf(dt)
    β = sup(dt)
    Xk = array(fp)
    for (i, X) in enumerate(Xk)
        if !isempty(tspan(X) ∩ dt)
            push!(idx, i)
        end
    end

    if isempty(idx)
        throw(ArgumentError("the time interval $dt does not intersect the time span, " *
                            "$(tspan(fp)), of the given flowpipe"))
    end
    return view(Xk, idx)
end

# concrete projection of a flowpipe along variables `vars`
function project(fp::Flowpipe, vars::NTuple{D, T}) where {D, T<:Integer}
    Xk = array(fp)
    πfp = map(X -> project(X, vars), Xk)
    return Flowpipe(πfp, fp.ext)
end

project(fp::Flowpipe, vars::Int) = project(fp, (vars,))
project(fp::Flowpipe, vars::AbstractVector{<:Int}) = project(fp, Tuple(vars))
project(fp::Flowpipe; vars) = project(fp, Tuple(vars))
project(fp::Flowpipe, i::Int, vars) = project(fp[i], vars)

# concrete linear map of a flowpipe for a given matrix
function linear_map(M::AbstractMatrix, fp::Flowpipe)
    Xk = array(fp)
    πfp = map(X -> linear_map(M, X), Xk)
    return Flowpipe(πfp, fp.ext)
end

"""
    shift(fp::Flowpipe{N, <:AbstractReachSet}, t0::Number) where {N}

Return the time-shifted flowpipe by the given number.

### Input

- `fp` -- flowpipe
- `t0` -- time shift

### Output

A new flowpipe such that the time-span of each constituent reach-set has been
shifted by `t0`.

### Notes

See also `Shift` for the lazy counterpart.
"""
function shift(fp::Flowpipe{N, <:AbstractReachSet}, t0::Number) where {N}
    return Flowpipe([shift(X, t0) for X in array(fp)], fp.ext)
end

"""
    convexify(fp::Flowpipe{N, <:AbstractLazyReachSet}) where {N}

Return a reach-set representing the convex hull array of the flowpipe.

### Input

- `fp` -- flowpipe

### Output

A reach-set that contains the convex hull array, `ConvexHullArray`, of the given
flowpipe.

### Notes

The time span of this reach-set is the same as the time-span of the flowpipe.

This function allocates an array to store the sets of the flowpipe.
"""
function convexify(fp::Flowpipe{N, <:AbstractLazyReachSet}) where {N}
    Y = ConvexHullArray([set(X) for X in array(fp)])
    return ReachSet(Y, tspan(fp))
end

"""
    convexify(fp::AbstractVector{<:AbstractLazyReachSet{N}}) where {N}

Return a reach-set representing the convex hull array of the array of the array of
reach-sets.

### Input

- `fp` -- array of reach-sets

### Output

A reach-set that contains the convex hull array, `ConvexHullArray`, of the given
flowpipe.

### Notes

The time span of this reach-set corresponds to the minimum (resp. maximum) of the
time span of each reach-set in `fp`.

This function allocates an array to store the sets of the flowpipe.

The function doesn't assume that the reach-sets are time ordered.
"""
function convexify(fp::AbstractVector{<:AbstractLazyReachSet{N}}) where {N}
    Y = ConvexHullArray([set(X) for X in fp])
    ti = minimum(tstart, fp)
    tf = maximum(tend, fp)
    return ReachSet(Y, TimeInterval(ti, tf))
end

# the dimension of sparse flowpipes is known in the type
function LazySets.dim(::Flowpipe{N, SparseReachSet{N, ST, D}}) where {N, ST, D}
    return D
end

tstart(F::Flowpipe, arr::UnitRange) = tstart(view(array(F), arr), contiguous=true)
tstart(F::Flowpipe, arr::AbstractVector) = tstart(view(array(F), arr))
tend(F::Flowpipe, arr::UnitRange) = tend(view(array(F), arr), contiguous=true)
tend(F::Flowpipe, arr::AbstractVector) = tend(view(array(F), arr))
tspan(F::Flowpipe, arr::UnitRange) = tspan(view(array(F), arr), contiguous=true)
tspan(F::Flowpipe, arr::AbstractVector) = tspan(view(array(F), arr))

# further setops
LazySets.is_intersection_empty(F::Flowpipe{N, <:AbstractLazyReachSet}, Y::LazySet) where {N} = all(X -> _is_intersection_empty(X, Y), array(F))
Base.:⊆(F::Flowpipe, X::LazySet) = all(R ⊆ X for R in F)
Base.:⊆(F::Flowpipe, Y::AbstractLazyReachSet) = all(R ⊆ set(Y) for R in F)

# lazy projection of a flowpipe
function Projection(F::Flowpipe, vars::NTuple{D, T}) where {D, T<:Integer}
    Xk = array(F)
    out = map(X -> Projection(X, vars), Xk)
    return Flowpipe(out)
end
Projection(F::Flowpipe; vars) = Projection(F, Tuple(vars))
Projection(F::Flowpipe, vars::AbstractVector{M}) where {M<:Integer} = Projection(F, Tuple(vars))

# membership test
function ∈(x::AbstractVector{N}, fp::Flowpipe{N, <:AbstractLazyReachSet{N}}) where {N}
    return any(R -> x ∈ set(R), array(fp))
end

function ∈(x::AbstractVector{N}, fp::VT) where {N, RT<:AbstractLazyReachSet{N}, VT<:AbstractVector{RT}}
    return any(R -> x ∈ set(R), fp)
end

function LazySets.linear_map(M, fp::Flowpipe)
    out = [linear_map(M, R) for R in fp]
    return Flowpipe(out, fp.ext)
end

function LazySets.affine_map(M, b, fp::Flowpipe)
    out = [reconstruct(R, affine_map(M, set(R), b)) for R in fp]
    return Flowpipe(out, fp.ext)
end

# lazy intersection
# use the flag `filter` (optional, default: `true`) to remove empty sets
@commutative function Base.:∩(fp::Flowpipe, X::LazySet; filter=true)
    out = [reconstruct(R, set(R) ∩ X) for R in fp]

    if filter
        # keep the sets which are not empty
        filter!(!isempty, out)
    end
    return Flowpipe(out, fp.ext)
end

# concrete intersection
# use the flag `filter` (optional, default: `true`) to remove empty sets
@commutative function LazySets.intersection(fp::Flowpipe, X::LazySet; filter=true)
    out = [reconstruct(R, intersection(set(R), X)) for R in fp]

    if filter
        # keep the sets which are not empty
        filter!(!isempty, out)
    end
    return Flowpipe(identity.(out), fp.ext)
end

# --------------------------------------------
# Specialized methods for template flowpipes
# --------------------------------------------

# assumes that the first reach-set is representative
function support_function_matrix(fp::Flowpipe{N, <:TemplateReachSet}) where {N}
    return support_function_matrix(first(fp))
end

# it is assumed that rows = (idx_pos_dir, idx_neg_dir) is such that each integer
# idx_pos_dir and idx_neg_dir refers to the row of the support function matrix
# correponding to a positive (resp. negative) direction
function flatten(fp::Flowpipe{N, <:TemplateReachSet}, rows=(1, 2)) where {N}
    # get the matrix of support function evaluations
    mat = support_function_matrix(fp)

    @assert size(mat, 1) ≥ 2 || throw(ArgumentError("the number of rows of the support function matrix should be at least 2, got $(size(mat, 1))"))
    @assert length(rows) == 2 || throw(ArgumentError("expected the number of rows of the support function matrix to be 2, got $(length(rows))"))
    idx_pos_dir = rows[1]
    idx_neg_dir = rows[2]

    RT = ReachSet{N, Interval{N, IA.Interval{N}}}
    out = Vector{RT}(undef, length(fp))

    @inbounds for (k, Rk) in enumerate(fp)
        Ik = Interval(-mat[idx_neg_dir, k], mat[idx_pos_dir, k])
        dtk = tspan(Rk)
        out[k] = ReachSet(Ik, dtk)
    end
    return Flowpipe(out, fp.ext)
end
