# ======================================
# Kronecker powers
# ======================================

using LinearAlgebra: checksquare

# TODO refactor to CarlemanLinearization.jl
import CarlemanLinearization: lift_vector

lift_vector(X0::IA.Interval, N) = lift_vector(Interval(X0), N)

"""
    kron_pow(x::IA.Interval, pow::Int)

Given an interval x and an integer pow, compute `x^pow`.

### Input

- `x`   -- interval
- `pow` -- integer

### Output

An interval enclosure of `x^pow`.

### Examples

```julia
julia> [kron_pow(2 .. 3, i) for i in 1:3]
3-element Array{IntervalArithmetic.Interval{Float64},1}:
  [2, 3]
  [4, 9]
 [8, 27]
```
"""
function kron_pow(x::IA.Interval, pow::Int)
    if pow == 1
        return x
    end

    return x^pow
end

"""
    kron_pow(x::Interval, pow::Int)

Given an interval x and an integer pow, compute `x^pow`.

### Input

- `x`   -- interval
- `pow` -- integer

### Output

An interval enclosure of `x^pow` as a LazySets `Interval`.
"""
function kron_pow(x::Interval, pow::Int)
    if pow == 1
        return x
    end

    return Interval(kron_pow(x.dat, pow))
end

"""
    kron_pow(H::AbstractHyperrectangle, pow::Int)

Given hyperrectangular set `H` and an integer `pow`, compute the Kronecker power
`H^{⊗ pow}`.

### Input

- `H`   -- hyperrectangular set
- `pow` -- integer power

### Output

A hyperrectangle.

### Algorithm

We compute `H^{⊗ pow}` where `H` is a hyperrectangular set by working with `H`
as a product of intervals.

See also `kron_pow` which requires `DynamicPolynomials.jl`.
"""
function kron_pow(H::AbstractHyperrectangle, pow::Int, algorithm=nothing)
    if pow == 1
        return H
    end

    if isnothing(algorithm)
        if isdefined(@__MODULE__, :DynamicPolynomials)
            algorithm = "symbolic"
        else
            algorithm = "explicit"
        end
    end
    if algorithm == "explicit"
        return _kron_pow_explicit(H, pow)
    elseif algorithm == "symbolic"
        return _kron_pow_symbolic(H, pow)
    else
        throw(ArgumentError("algorithm \"$algorithm\" is not known"))
    end
end

function _kron_pow_explicit(H::AbstractHyperrectangle, pow::Int)
    x = [Interval(low(H, j), high(H, j)) for j in 1:dim(H)]
    r = reduce(kron, fill(x, pow))

    low_r = min.(r)
    high_r = max.(r)
    return Hyperrectangle(; low=low_r, high=high_r)
end

"""
    kron_pow_stack(x::IA.Interval, pow::Int)

Return a hyperrectangle with the interval powers `[x, x^2, …, x^pow]`.

### Input

- `x`   -- interval
- `pow` -- integer power

### Output

A hyperrectangle such that the `i`-th dimension is the interval `x^i`.
"""
function kron_pow_stack(x::IA.Interval, pow::Int)
    return convert(Hyperrectangle, IntervalBox([kron_pow(x, i) for i in 1:pow]))
end

"""
    kron_pow_stack(x::Interval, pow::Int)

Return a hyperrectangle with the interval powers `[x, x^2, …, x^pow]`.

### Input

- `x`   -- interval
- `pow` -- integer power

### Output

A hyperrectangle such that the `i`-th dimension is the interval `x^i`.
"""
function kron_pow_stack(x::Interval, pow::Int)
    return convert(Hyperrectangle, CartesianProductArray([kron_pow(x, i) for i in 1:pow]))
end

"""
    kron_pow_stack(H::AbstractHyperrectangle, pow::Int)

Return the Cartesian product array ``H × H^{⊗2} × ⋯ × H^{⊗pow}`` where ``H`` is a
hyperrectangular set and ``H^{⊗ i}`` is the `i`-th Kronecker power of ``H``.

### Input

- `H`   -- hyperrectangular set
- `pow` -- integer power

### Output

A Cartesian product array of hyperrectangles.
"""
function kron_pow_stack(H::AbstractHyperrectangle, pow::Int)
    out = [kron_pow(H, i) for i in 1:pow]
    return CartesianProductArray(out)
end

# ------------------------------------------------------------
# Functionality that requires DynamicPolynomials.jl
# ------------------------------------------------------------

function load_kron_dynamicpolynomials()
    return quote
        function _kron_pow_symbolic(H::AbstractHyperrectangle{N}, pow::Int) where {N}
            n = dim(H)
            @polyvar x[1:n]
            B = convert(IntervalBox, H)
            dict = Dict((x[i] => i for i in 1:n)...)
            y = kron_pow(x, pow)

            out = Vector{IA.Interval{N}}(undef, length(y))
            for (i, p) in enumerate(y)
                aux = interval(1)
                for (xj, j) in powers(p)
                    aux = aux * B[dict[xj]]^j
                end
                out[i] = aux
            end
            Bpow = IntervalBox(out)
            Hpow = convert(Hyperrectangle, Bpow)
            return Hpow
        end
    end
end  # quote / load_kron_dynamicpolynomials()

# ------------------------------------------------------------
# Functionality that requires MultivariatePolynomials.jl
# ------------------------------------------------------------

function load_kron_multivariate()
    return quote
        """
            kron_pow(x::Vector{<:AbstractVariable}, pow::Int)

        Compute the higher order concrete Kronecker power: `x ⊗ x ⊗ ... ⊗ x`, `pow` times
        for a vector of symbolic monomials.

        ### Input

        - `x`   -- polynomial variable
        - `pow` -- integer

        ### Output

        Vector of multivariate monomial corresponding to `x^{⊗ pow}`.

        ### Examples

        ```julia
        julia> using DynamicPolynomials

        julia> @polyvar x[1:2]
        (PolyVar{true}[x₁, x₂],)

        julia> x
        2-element Array{PolyVar{true},1}:
         x₁
         x₂

        julia> kron_pow(x, 2)
        4-element Array{Monomial{true},1}:
         x₁²
         x₁x₂
         x₁x₂
         x₂²
        ```
        """
        function kron_pow(x::Vector{<:AbstractVariable}, pow::Int)
            @assert pow > 0 "expected positive power, got $pow"
            if pow == 1
                return x
            else
                return kron(x, kron_pow(x, pow - 1))
            end
        end
    end
end  # quote / load_kron_multivariate()
