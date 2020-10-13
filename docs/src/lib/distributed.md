DocTestSetup = :(using ReachabilityAnalysis)
```

# Distributed computations

```@contents
Pages = ["distributed.md"]
Depth = 3
```

This section of the manual describes functions to make use of distributed computation.

```@contents
Pages = ["distributed.md"]
```

```@meta
CurrentModule = ReachabilityAnalysis
```

## Parallel solve

Reachability computations can be computed in parallel using Julia's built-in
multithreaded parallelism. This feature is available formulating and solving
an initial-value problem with an array of initial conditions.

`

## About BLAS threads

To control the number of threads used by your BLAS library, use the function
`Base.LinAlg.BLAS.set_num_threads(n)`, where `n` is an integer. Furthermore,
the function `get_num_threads()` defined below will return the current value.

*Note.* If you are using Julia v"0.7-" (run the command `VERSION` to find this),
instead of `Base.LinAlg` below use `LinearAlgebra`, and this module should have
been loaded in the current scope with `using LinearAlgebra`.

```julia
#
# This function is a part of Julia. License is MIT: https://julialang.org/license
#
function get_num_threads() # anonymous so it will be serialized when called
    blas = Base.LinAlg.BLAS.vendor()
    # Wrap in a try to catch unsupported blas versions
    try
        if blas == :openblas
            return ccall((:openblas_get_num_threads, Base.libblas_name), Cint, ())
        elseif blas == :openblas64
            return ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
        elseif blas == :mkl
            return ccall((:MKL_Get_Max_Num_Threads, Base.libblas_name), Cint, ())
        end

        # OSX BLAS looks at an environment variable
        if Sys.isapple()
            return ENV["VECLIB_MAXIMUM_THREADS"]
        end
    end

    return nothing
end
```
