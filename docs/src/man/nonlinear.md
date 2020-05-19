# Nonlinear ordinary differential equations

In this section we illustrate the flowpipe computation for a nonlinear system.

## Model description

Our running example is the [Lotka-Volterra model](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations).
The 2-dimensional Lotka-Volterra system depicts the populations change of a class of predators and a class of
preys. The growth rate of preys’ population $x$ over time is given by

```math
\dot{x} = x\cdot (\alpha - \beta \cdot y)
```
wherein  $\alpha, \beta$ are constant parameters and $y$ is the population of predators.

It gives that the number of preys grows exponentially without predation.

The population growth of predators is governed by the differential equation

```math
\dot{y} = -y\cdot (\gamma - \delta\cdot x)
```
wherein  $\gamma, \delta$ are constant parameters.

We set those parameters as  $\alpha = 1.5 ,  \beta = 1 ,  \gamma = 3$  and  $\delta = 1$.

```@example lotka_volterra
using ReachabilityAnalysis

@taylorize function lotka_volterra!(du, u, p, t)
    local α, β, γ, δ = 1.5, 1.0, 3.0, 1.0
    du[1] = u[1] * (α - β*u[2])
    du[2] = -u[2] * (γ - δ*u[1])
    return du
end
```

## Reachability settings

The reachability settings are taken from [this resource](https://ths.rwth-aachen.de/research/projects/hypro/lotka-volterra/).

We consider the initial set  $x\in [4.8,5.2], y \in [1.8,2.2]$.

```@example lotka_volterra
X₀ = Hyperrectangle(low=[4.8, 1.8], high=[5.2, 2.2])

prob = @ivp(x' = lotka_volterra!(x), dim: 2, x(0) ∈ X₀)
```

## Results

We compute the flowpipe using the TMJets algorithm for the time horizon $[0,5]$:

```@example lotka_volterra
sol = solve(prob, T=5.0)

setrep(sol)
```

We can change to the zonotopic overapproximation of the flowpipe using
the `overapproximate` function:

```@example lotka_volterra
sol = overapproximate(sol, Zonotope)

setrep(sol)
```

Finally we plot the solution in phase-space:

```@example lotka_volterra
using Plots

plot(sol, vars=(1, 2), xlab="x", ylab="y", lw=0.2, color=:lightblue, lab="Flowpipe")
plot!(X₀, color=:orange, lab="Xo")
```

## Some common gotchas

### What is `@taylorize`? Do I need it?

`@taylorize` is a macro which parses the functions containing the ODEs to be integrated,
allowing to speed up repeated evaluations. The macro is defined in
[TaylorIntegration.jl](https://github.com/PerezHz/TaylorIntegration.jl), see
[`@taylorize`'s documentation in TaylorIntegration.jl](https://perezhz.github.io/TaylorIntegration.jl/stable/taylorize/)
for further details. Since it is an optimization, it is *not* mandatory,
though it is recommended as it helps to reduce the number of allocations and as a
consequence it usually gives a performance boost.

### How can I get the most of out `@taylorize`?

The main advice is to refactor expressions involving several terms into smaller
expressions which involve at most two arguments, making appropriate use of parentheses
if needed. For further limitations and advice see [this section of TaylorInegrations.jl's
documentation](https://perezhz.github.io/TaylorIntegration.jl/stable/taylorize/#Limitations-and-some-advices-1).

**Example.** Here is an example that uses some of the above recommendations.
Start with `f!` defined below:

```julia
@taylorize function f!(du, u, params, t)
    local a = 0.3
    x, y, z = u[1], u[2], u[3]

    du[1] = -x * y/(1 + x)
    du[2] = x * y/(1 + x) - a * y
    du[3] = a * y * y
    return du
end
```

Observe that the terms `x * y` can be factored out into a new auxiliary variable
`aux`, and all intermediate operations involve two arguments at most using parentheses:

```julia
@taylorize function g!(du, u, params, t)
    local a = 0.3
    x, y, z = u[1], u[2], u[3]

    num = x * y
    den = 1 + x
    aux = num/den
    du[1] = -aux
    du[2] = aux - a * (y * y)
    du[3] = a * (y * y)
    return du
end
```

### Why do I get an `ArgumentError` when trying to plot a solution?

For flowpipes computed using algorithm `TMJets`, the representation is a
*Taylor model flowpipe*, and as such we don't have methods to visualize them exactly.
However, to reason about such flowpipes, either for plotting or performing set-based
operations, we can *overapproximate* them with other set representations
-- usually, convex sets such as boxes or zonotopes--. For instance, the command
`overapproximate(sol, Zonotope)` applies `overapproximate(Ri, Zonotope)` for each
reach-set `Ri` in the solution `sol`.

!!! note
    The default plotting behavior may change in the future, see discussion in
    issue [#173](https://github.com/JuliaReach/ReachabilityAnalysis.jl/issues/173).

**Example.** Consider the differential equation `x'(t) = x^2(t) - 1`, given the interval
initial condition `[0.4, 0.5]`:

```julia
function f!(dx, x, p, t)
    dx[1] = x[1]^2 - 1.0
end
prob = @ivp(x' = f!(x), dim: 1, x(0) ∈ 0.4 .. 0.5)
sol = solve(prob, T=10.0)
```
Trying to plot the solution with the command `plot(sol, vars=(1, 2))` will fail with
an `ArgumentError`. However, you can plot the zonotopic overapproximation of this flowpipe:

```julia
solz = overapproximate(sol, Zonotope)
plot(solz, vars=(0, 1))
```

### Equations with constant terms (`BoundsError`)

Equations that involve constant terms may give a `BoundsError`. This is a known bug
(cf. issue [#179](https://github.com/JuliaReach/ReachabilityAnalysis.jl/issues/179))
and it is related to cases in which the update rule for the right-hand side does not
have the expected coefficient type. The current solution is to change terms
like `du[1] = 1.0` to `du[1] = 1.0 + zero(u[1])`, i.e. let Julia's promotion mechanism
take care by adding the given numeric constant with the zero element of the type of `u`.

**Example.** Consider the function `f!`:

```julia
@taylorize function f!(du, u, p, t)
    du[1] = u[3]^3 - u[2] + u[4]
    du[2] = u[3]
    du[3] = 2.0
    du[4] = u[4]
    return du
end
```
Integrating this function will likely fail with a `BoundsError`. However, we can
re-write it in this way:

```julia
@taylorize function f!(du, u, p, t)
    local two = 2.0 + zero(u[1])
    du[1] = u[3]^3 - u[2] + u[4]
    du[2] = u[3]
    du[3] = two
    du[4] = u[4]
    return du
end
```
