"""
    AbstractPost

Abstract supertype of all post operator types.

### Notes

All post operators should provide the following methods:

```julia
init(op::AbstractPostOperator, system, options)

post(op::AbstractPostOperator, system, ...)
```
"""
abstract type AbstractPost end

"""
    AbstractContinuousPost

Abstract supertype of all continuous post operators.
"""
abstract type AbstractContinuousPost <: AbstractPost end

"""
    AbstractDiscretePost

Abstract supertype of all discrete post operators.

### Notes

All discrete post operators should provide the following method, in addition
to those provided for general post operators:

```julia
tube⋂inv!(𝒫::AbstractDiscretePost, reach_tube::Vector{<:AbstractReachSet{<:LazySet, N}},
          invariant, Rsets, start_interval) where {N}
```
"""
abstract type AbstractDiscretePost <: AbstractPost end
