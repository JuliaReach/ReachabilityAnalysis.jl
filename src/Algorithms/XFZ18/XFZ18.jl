"""
References:

https://dl.acm.org/doi/10.1145/3178126.3178133
"""
struct XFZ18{N, ST} <: AbstractContinuousPost
#
end

numtype(::XFZ18{N}) where {N} = N
