```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Solutions

## Abstract interface

```@docs
ReachabilityAnalysis.AbstractSolution
```

## Solution of a reachability problem

```@docs
ReachabilityAnalysis.ReachSolution
```

## Solution of a verification problem

```@docs
ReachabilityAnalysis.CheckSolution
```

## Solving a reachability problem

```@docs
solve
AbstractPost
AbstractContinuousPost
AbstractDiscretePost
```

## Solving a hybrid reachability problem

```@docs
AbstractWaitingList
WaitingList
StateInLocation
```

TODO: document other methods in `solutions.jl`.
