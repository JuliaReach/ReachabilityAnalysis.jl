"""
    ORBIT{N, VT, AM} <: AbstractContinuousPost

Implementation of discrete-time integration for deterministic linear ODEs
for singleton initial conditions and singleton input sets.

## Fields

- `δ`            -- step-size of the discretization
- `approx_model` -- (optional, default: `NoBloating`) approximation model

The type fields are:

- `N`   -- number type of the step-size
- `AM`  -- type of the approximation model

## Notes

Reachability solutions computed with the `ORBIT` algorithm have its own plot
recipe, producing a scatter plot. Use the `markershape` option to change the shape
of the markers used for the scatter plot. To "connect" the points, use the option
`seriestype=:path`. For additional options see the [Plots.jl documentation](https://docs.juliaplots.org/latest/generated/attributes_series/).
"""
struct ORBIT{N, VT, AM} <: AbstractContinuousPost
    δ::N
    approx_model::AM
end

step_size(alg::ORBIT) = alg.δ
numtype(::ORBIT{N}) where {N} = N

# convenience constructor using symbols
function ORBIT(; δ::N, approx_model::AM=NoBloating(exp=:base, setops=:concrete)) where {N, AM}
    VT = Vector{N}
    return ORBIT{N, VT, AM}(δ, approx_model)
end

function setrep(::ORBIT{N, VT, AM}) where {N, VT, AM}
    Singleton{N, VT}
end

function rsetrep(alg::ORBIT{N}) where {N}
    ST = sertrep(alg)
    ReachSet{N, ST}
end

include("post.jl")
