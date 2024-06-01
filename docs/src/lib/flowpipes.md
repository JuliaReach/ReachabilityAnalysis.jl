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
LazySets.ρ(::AbstractVector, ::AbstractFlowpipe)
LazySets.σ(::AbstractVector, ::AbstractFlowpipe)
LazySets.dim(::AbstractFlowpipe)
set(::AbstractFlowpipe, ::Integer)
set(::AbstractFlowpipe, ::AbstractVector)
set(::AbstractFlowpipe)
tstart(::AbstractFlowpipe)
tend(::AbstractFlowpipe)
tspan(::AbstractFlowpipe)
vars(::AbstractFlowpipe)
```

## Flowpipe

```@docs
Flowpipe
location(::Flowpipe)
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

## Mapped flowpipe

```@docs
MappedFlowpipe
```
