```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

## Benchmark repository

The benchmark suite of JuliaReach is available at the
[ReachabilityBenchmarks](https://github.com/JuliaReach/ReachabilityBenchmarks) repository.

The repository includes:

- Benchmarks for all models in the `ReachabilityAnalysis.jl` test suite.
- Benchmarks for models which are not part of the `ReachabilityAnalysis.jl` test suite,
  and are used to detect regressions (i.e. examples that accidentally run slower,
  due to changes in core downstream libraries such as `LazySets.jl`).
- The [SLICOT models](http://slicot.org/matlab-toolboxes/model-reduction), from the
  SLICOT Model and Controller Reduction Toolbox, which reflect real world applications
  with dimensions ranging from several dozens to over 10.000.
