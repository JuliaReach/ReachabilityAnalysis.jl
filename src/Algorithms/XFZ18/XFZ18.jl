"""
Implementation of the reachability algorithm for the class of polynomial ODEs
with uncertain initial states (see [XueFZ18](@citet)).
This method consists of reducing the Hamilton-Jacobi-Bellman equation to a
hierarchy of semidefinite programs that are solved using an SDP solver.

### Algorithm

We refer to [XueFZ18](@citet) for technical details.
"""
struct XFZ18{N,ST} <: AbstractContinuousPost
    #
end

numtype(::XFZ18{N}) where {N} = N
