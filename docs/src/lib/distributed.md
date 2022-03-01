```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Distributed computations

This section of the manual gives additional pointers on making use of distributed computation
within the library.

```@contents
Pages = ["distributed.md"]
Depth = 3
```

## Multi-threading

Since version 1.3, Julia has built-in multi-threaded parallelism and we have
applied such feature at different levels in the library:

1. Multi-threaded implementations of reachability algorithms. This approach
   relies on the particular structure of each algorithm and hence it doesn't
   generalize to different algorithms. An example is the exploration
   of template directions in simultaneous with `LGG09`. Details are included in the
   documentation page of each algorithm.

2. Parallel `solve`, with split initial conditions. This approach applies to all
   algorithms, because we can always split an initial domain into smaller subdomains,
   and run `solve` on each of these regions in parallel. See [Parallel solve](@ref) below
   for details.

[Julia's documentation on multi-threading](https://docs.julialang.org/en/v1/manual/multi-threading/)
describes how to check the number of threads available in the current session, `Threads.nthreads()`,
and how to control the number of threads e.g. by starting Julia with `$ julia --threads 4` to
use four threads. Please note that to make the number of threads persistent across different
Julia sessions you should export an environment variable, `export JULIA_NUM_THREADS=4` to be
added in your `.bashrc` or `.bash_profile` files, or setting `ENV["JULIA_NUM_THREADS"]=4`
in your `.julia/config/startup.jl` file.

## GPGPU computing

The recommended entry point for using general-purpose graphical processing units (GPGPU)
in Julia is the library [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl).
The package exports an array type `CuArray` used as an abstraction to perform array
operations on the GPU device; common linear algebra operations are readily available
through dispatch on `CuArray`. Apart from convenient high-level syntax to do
basic linear algebra on GPUs; more advanced kernel functions can be implemented as well
as it is explained in the [CUDA.jl documentation](https://juliagpu.gitlab.io/CUDA.jl/).

## Parallel solve

Reachability computations can be computed in parallel using Julia's built-in
multithreaded parallelism. This feature is available formulating and solving
an initial-value problem with an array of initial conditions.

## About BLAS threads

To control the number of threads used by your BLAS library, use the function
`BLAS.set_num_threads(n)`, where `n` is an integer. Furthermore,
the function `BLAS.get_num_threads()` defined below will return the current value.

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
