"""
    BOX{N, AM, D} <: AbstractContinuousPost

Implementation of a reachability method for linear systems using box approximations.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default: `Forward`) approximation model;
                    see `Notes` below for possible options
- `dim`          -- (optional default: `missing`) ambient dimension

## Notes

The type fields are:

- `N`         -- number type of the step-size
- `AM`        -- approximation model
- `D`         -- refers to the dimension of the system
- `recursive` -- (optional default: `false`) if `true`, use the implementation that
                 recursively computes each reach-set; otherwise, use the implementation
                 that unwraps the sequence until the initial set

The default approximation model is:

```julia
Forward(sih=:concrete, exp=:base, phi2=:base, setops=:lazy)
```

This algorithm solves the set-based recurrence equation ``X_{k+1} = ΦX_k ⊕ V_k``
by computing a tight hyperrectangular over-approximation of ``X_{k+1}`` at each
step ``k ∈ \\mathbb{N}``. The recursive implementation uses the previously computed
set ``X_k`` to compute ``X_{k+1}``. However, it is known that this method incurs
wrapping effects. The non-recursive implementation instead computes ``X_{k+1}``
by unwrapping the discrete recurrence until ``X_0 = Ω₀``, at the expense of computing
powers of the matrix ``Φ``. These ideas are discussed in [[BFFPSV18]](@ref).

### References

This algorithm is essentially a non-decomposed version of the method in [[BFFPSV18]](@ref),
using hyperrectangles as set representation. For a general introduction we refer
to the dissertation [[LG09]](@ref).

Regarding the approximation model, by default we use an adaptation of the method
presented in [[FRE11]](@ref).
"""
@with_kw struct BOX{N, AM, D} <: AbstractContinuousPost
    δ::N
    approx_model::AM=Forward(sih=:concrete, exp=:base, phi2=:base, setops=:lazy)
    static::Bool=false
    dim::D=missing
    recursive::Bool=false
end

step_size(alg::BOX) = alg.δ
numtype(::BOX{N}) where {N} = N

function rsetrep(alg::BOX{N}) where {N}
    if !alg.static
        VT = Vector{N}
    else
        n = alg.dim
        VT = SVector{n, N}
    end
    RT = ReachSet{N, Hyperrectangle{N, VT, VT}}
    return RT
end

include("post.jl")
include("reach_homog.jl")
include("reach_inhomog.jl")
