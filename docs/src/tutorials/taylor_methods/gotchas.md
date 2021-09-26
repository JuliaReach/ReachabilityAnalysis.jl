```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Some common gotchas

In this section we comment on some performance aspects about the
solution of initial-value problems for nonlinear systems using Taylor methods.

### What is `@taylorize`? Do I need it?

`@taylorize` is a macro which parses the functions containing the ODEs to be integrated,
allowing to speed up repeated evaluations. The macro is defined in
[TaylorIntegration.jl](https://github.com/PerezHz/TaylorIntegration.jl), see
[`@taylorize`'s documentation in TaylorIntegration.jl](https://perezhz.github.io/TaylorIntegration.jl/stable/taylorize/)
for further details. Since it is an optimization, it is *not* mandatory,
though it is recommended as it helps to reduce the number of allocations and as a
consequence it usually gives a performance boost.

### How can I get the most out of `@taylorize`?

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
`aux`, and all intermediate operations can be arranged to only involve two arguments,
using parentheses:

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

### How are solutions obtained with Taylor models methods plotted?

Flowpipes computed using algorithm `TMJets` (or its variations),
use Taylor model reach-set representations (`TaylorModelReachSet`),
which define an implicit set in time and in space. Since exact visualization of
such objects is difficult (and often unnecessary), we resort to an outer approximation
with simpler sets. Either for plotting or performing set-based operations, we
can *overapproximate* a `TaylorModelReachSet` with other set representations -- usually,
convex sets such as boxes or zonotopes--. The command `overapproximate(sol, Zonotope)`
applies `overapproximate(Ri, Zonotope)` for each reach-set `Ri` in the solution `sol`.

By default, when plotting the solution obtained with such solvers, the zonotopic
overapproximation of the flowpipe is used, with a single zonotope
per Taylor model reach-set. Such approximation, while it is generally coarse,
is often sufficient for visualization purposes.

### Equations with constant terms (`BoundsError`)

Equations that involve constant terms may give a `BoundsError`. This is a known bug
(cf. issue [#179](https://github.com/JuliaReach/ReachabilityAnalysis.jl/issues/179))
and it is related to cases in which the update rule for the right-hand side does not
have the expected coefficient type. The current solution is to change terms
like `du[1] = 1.0` into `du[1] = 1.0 + zero(u[1])`, i.e. let Julia's promotion mechanism
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
