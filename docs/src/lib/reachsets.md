```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Reach-sets

## Abstract interface

```@docs
AbstractReachSet
```

The functions are available at the interface level.

```@docs
set(::AbstractReachSet)
setrep(::AbstractReachSet)
tspan(::AbstractReachSet)
tstart(::AbstractReachSet)
tend(::AbstractReachSet)
dim(::AbstractReachSet)
copy(::AbstractReachSet)
shift(::AbstractReachSet)
```

```@docs
AbstractLazyReachSet
project(::AbstractLazyReachSet, ::NTuple{D,Int}; check_vars::Bool=true) where {D}
```

## Reachable set

```@docs
ReachSet
```

## Sparse reachable set

```@docs
SparseReachSet
```

## Taylor model reach-sets

```@docs
AbstractTaylorModelReachSet
TaylorModelReachSet
```


## Template reach-sets

```@docs
TemplateReachSet
```
