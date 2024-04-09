```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Linear reachability methods

This section begins by introducing the notion of reachable set (*reach-set*).
We show how to visualize and perform set operations with them. Then we introduce *flowpipes*
as the union of reach-sets and illustrate with examples how to perform different
visualization and algebraic operations with flowpipes. Then we consider a simple
set-propagation problem that consists of a two-dimensional rotation.

After those preliminary sections, we discuss the notion of conservative time
discretization for systems of linear differential equations of the form

```math
    x'(t) = Ax(t),\qquad x(0) ∈ X_0 ∈ \mathbb{R}^n.
```
for all times ``t ∈ [0, T]``. Linear systems with non-deterministic inputs
is discussed in another section of this manual. We also consider an invariant
specification: ``x(t) ∈ X`` for all times.

The final part of this section introduces support function techniques and discusses
the helicopter model application. For ease of exposition, this section only considers.
We also show how to introduce state-space invariants, i.e. to impose

## Quickstart guide

The user interface to solve initial-value problems is quite simple.

```@example quickstart
using ReachabilityAnalysis

# initial-value problem specification
p = @ivp(x' = -x, x(0) ∈ Interval(1, 2))

# flowpipe computation
sol = solve(p, T=5)

nothing # hide
```

Here we have solved the differential equation

```math
x'(t) = -x(t),\qquad x(0) ∈ X_0 = [1, 2] ⊆ \mathbb{R},
```
for ``t ∈ [0, 5]``, whose solution is known to be the decaying exponential
``x(t) = x_0e^{-t}``. Let's plot the result, which requires loading the [Plots.jl]()
package

variable index `0` denotes time).

```@example quickstart
using Plots

# post-processing or plotting
plot(sol, vars=(0, 1), xlab="t", ylab="x(t)")
```

For comparison we can plot the trajectories at the endpoints of
the initial interval. The example also shows that the black line segments
can be hidden by setting the `lw=0` argument (for `linewidth`).

```@example quickstart
trange = range(0, 5, length=100)

plot(sol, vars=(0, 1), xlab="t", ylab="x(t)", lw=0)
plot!(trange, 1.0 * exp.(-trange), vars=(0, 1), c=:magenta, lab="")
plot!(trange, 2.0 * exp.(-trange), vars=(0, 1), c=:magenta, lab="")
```

The user-facing interface is designed to be intuitive and interactive, and it is
inspired by of other Julia packages such as [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).
Moreover, the library internals are written in a modular and composable way,
such that advanced users are able to modify and changed easily, or to compose with
other algorithms, different steps of the solution process.

In the rest of this notebook we explore different problem specifications,
different algorithm choices as well as some processing capabilities. As a preliminary,
we introduce the concept of reach-set in the context of the toy model presented above.

## Introducing reach-sets

Basically, a *reach-set* is a structure that holds a set and a time span associated to it.
Here we pick the first reach-set and show that it is of an interval type.

```@example quickstart
# solutions implement the array interface
R = sol[1]
```
Note that `sol[1]` works since solution structures implement Julia's array interface
-- meaning that slicing also works, and e.g. `sol[end-3:end]` returns the last
three reach-sets computed.

The solution obtained by set propagation consists of a *flowpipe*, which is just an
array of reach-sets, and behaves like their set union. Flowpipes are at the right level
of abstraction concerning solutions obtained with set propagation methods.

We can plot the first reach-set as well.

```@example quickstart
trange = range(0, 5, length=100)

plot(sol, vars=(0, 1), xlab="t", ylab="x(t)", lw=0)
plot!(R, vars=(0, 1), xlab="t", ylab="x(t)", lw=0, alpha=1., c=:orange, lab="R = sol[1]")
plot!(trange, 1.0 * exp.(-trange), vars=(0, 1), c=:magenta, lab="")
plot!(trange, 2.0 * exp.(-trange), vars=(0, 1), c=:magenta, lab="")
```

We can also plot the tenth reach-set in red:

```@example quickstart
plot!(sol[10], vars=(0, 1), xlab="t", ylab="x(t)", lw=0, alpha=1., c=:red, lab="sol[10]")
```

The type of reach-sets specifies its numeric type as well as the *set representation* used;
in this case, an interval.

```@example quickstart
typeof(R)
```
Observing the horizontal axis of the plots reveals that reach-sets have a non-zero
along on the time axis. By construction, it is guaranteed that the flowpipe is an
*enclosure* of the true solutions, i.e. no trajectory escapes each reach-set for a given
time-span. The time span associated to a reach-set is obtained with the function `tspan`.

```@example quickstart
tspan(R)
```

Actually, the algorithm that has been used to solve the problem has a fixed step size
of ``δ = 0.05``.

```@example quickstart
tspan(sol[end]) # time-span of the final reach-set
```

We can check by extracting the algorithm information from the solution struct:

```@example quickstart
# methods' step-size
sol.alg.δ
```

The functions `tstart` and `tend` return the starting and final time of the reach-set.

```@example quickstart
tstart(R)
```

```@example quickstart
tend(R)
```

The set wrapped by ``R`` is obtained with `set`:

```@example quickstart
set(R)
```

We see it is an interval (`LazySets.Interval` is just a thin wrapper around `IntervalArithmetic.Interval`;
more on this in a note below).

It is interesting to observe that the infimum of `sol[1]` 0.9487, while we have
specified that the initial interval is ``X_0 = [1, 2]``. The explanation is that the
computed reach-set contains the true solution for *all* intermediate times between
0 and 0.05 and for all initial states in ``X_0``. Since the values of ``x(t)`` decrease
in the time interval ``[0, 0.05]`` the width of `sol[1]` should be sufficiently big
as to enclose those variations. We can make a quick check:

```@example quickstart
1.0 * exp(-0.05)
```
which shows that `R` indeed is a correct enclosure of the solution for all points
in the initial set. It is easy to check that by decreasing the step-size, the enclosure
of the solution at ``[0, 0.05]`` converges to the true value (from outside), and similarly
for other time intervals.

### More about time intervals

The time span associated to a reach-sets is an *interval* satisfying the rules of interval arithmetic. We use `IntervalArithmetic.jl` to represent time intervals.

```@example quickstart
ReachabilityAnalysis.TimeInterval # alias
```

```@example quickstart
# retrieve the time span associated to this reach-set
tspan(R)
```

```@example quickstart
# it is an interval (from IntervalArithmetic.jl)
typeof(tspan(R))

```
Time intervals are represented using intervals from `IntervalArithmetic.jl`. This choice guarantees that all calculations involving time are carried out using rigorous floating-point calculations with interval arithmetic: all quantities are treated as intervals, which are propagated throughout a calculation. The final result is an interval that is guaranteed to contain the correct result, starting from the given initial data.

If desired, it is *also* work with set using interval elements:

```@example quickstart
Bint = BallInf(interval.(ones(2)), interval(0.2))
```
Finally, note that if you create a reach-set by passing a time point, it is automatically converted to an interval:

```@example quickstart
R = ReachSet(rand(BallInf, dim=5), 1.0)

tspan(R)
```
