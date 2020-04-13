```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Flowpipes

## Abstract interface

```@docs
AbstractFlowpipe
```

The following functions are available at the interface level.

```@docs
basetype(::Type{<:AbstractFlowpipe})
LazySets.ρ(::AbstractVector, ::AbstractFlowpipe)
LazySets.σ(::AbstractVector, ::AbstractFlowpipe)
LazySets.dim(::AbstractFlowpipe)
set(::AbstractFlowpipe, ::Integer)
tstart(::AbstractFlowpipe)
tend(::AbstractFlowpipe)
tspan(::AbstractFlowpipe)
```

## Flowpipe

```@docs
Flowpipe
```

The following methods are available.

```@docs
shift
```

## Hybrid flowpipe

```@docs
HybridFlowpipe
```


## Mixed flowpipe

```@docs
MixedFlowpipe
```

## Shifted flowpipe

```@docs
ShiftedFlowpipe
```


## Partitioned flowpipe

```@docs
PartitionedFlowpipe
```
