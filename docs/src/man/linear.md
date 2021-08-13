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
    x'(t) = Ax(t),\qquad x(0) \in X_0 \in \mathbb{R}^n.
```
for all times $t \in [0, T]$. Linear systems with non-deterministic inputs
is disussed in another section of this manual. We also consider an invariant specification: $x(t) \in X$ for all times.

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
x'(t) = -x(t),\qquad x(0) \in X_0 = [1, 2] \subset \mathbb{R},
```
for $t \in [0, 5]$, whose solution is known to be the decaying exponential
$x(t) = x_0e^{-t}$. Let's plot the result, which requires loading the [Plots.jl]()
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
inspired by of other Julia packages such as DifferentialEquations.jl.
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
of $δ = 0.05$.

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

The set wrapped by $R$ is obtained with `set`:

```@example quickstart
set(R)
```

We see it is an interval (`LazySets.Interval` is just a thin wrapper around `IntervalArithmetic.Interval`;
more on this in a note below).

It is interesting to observe that the infimum of `sol[1]` 0.9487, while we have
specified that the initial interval is $X_0 = [1, 2]$. The explanation is that the
computed reach-set contains the true solution for *all* intermediate times between
0 and 0.05 and for all initial states in $X_0$. Since the values of $x(t)$ decrease
in the time interval $[0, 0.05]$ the width of `sol[1]` should be sufficiently big
as to enclose those variations. We can make a quick check:

```@example quickstart
1.0 * exp(-0.05)
```
which shows that `R` indeed is a correct enclosure of the solution for all points
in the initial set. It is easy to check that by decreasing the step-size, the enclosure
of the solution at $[0, 0.05]$ converges to the true value (from outside), and similarly
for other time intervals.

### On time intervals representation

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

## Set propagation in discrete time

Our motivating example is to solve the following simple linear set-based recurrence

```math
X_{k+1} = M(\theta) X_k, \qquad 0 \leq k \leq 50, X_0 = [0.8, 1.2] \times [0.8, 1.2] \subset \mathbb{R}^2
```

Let $\theta \in [0, 2 \pi]$ be an equally spaced vector of length $50$, and $M(\theta)$ is the [rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix) given by:

```math
M(\theta) = \begin{pmatrix}\cos\theta && \sin\theta \\ -\sin\theta && \cos\theta \end{pmatrix}.
```

The matrix $M(\theta)$ rotates points in the xy-plane clockwise through an angle $θ$ around the origin of a two-dimensional Cartesian coordinate system.

### Propagating point clouds

To gain some intuition let's build the matrix and apply it to some points.

```@example discrete_propagation
using ReachabilityAnalysis # hide
using ReachabilityAnalysis: center  # hide

import Plots: plot, plot!, xlims!, ylims! # hide

# initial set
X0 = BallInf(ones(2), 0.2)

# rotation matrix
M(θ) = [cos(θ) sin(θ); -sin(θ) cos(θ)]
```

```@example discrete_propagation
plot(X0, c=:white)
plot!(M(pi/4) * X0, c=:white)

# center of the initial set
c = center(X0) |> Singleton

# list of vertices
V = vertices_list(X0) .|> Singleton |> UnionSetArray

plot!(c)
plot!(V)
plot!(M(pi/4) * c) # rotate 45 degrees
plot!(M(pi/4) * V)

xlims!(0.0, 1.8) # hide
ylims!(-0.4, 1.4) # hide
```

```@example discrete_propagation
X = sample(X0, 25) .|> Singleton |> UnionSetArray

plot!(X)
plot!(M(pi/4) * X, c=:blue)     # rotate 45 degrees
plot!(M(pi/2) * X, c=:green)    # rotate 90 degrees
plot!(M(4pi/3) * X, c=:orange)  # rotate 180 degrees

plot!(X0, c=:white)
plot!(M(pi/4) * X0, c=:white)
plot!(M(pi/2) * X0, c=:white)
plot!(M(4pi/3) * X0, c=:white)

xlims!(0.0, 1.8) # hide
ylims!(-0.4, 1.4) # hide
```

Does propagating point clouds solve the problem?
Practically speaking, while we can compute the successors of any
$x_0 \in X_0$ we still lack a *global* description of the set according
the the given discrete recurrence.. which brings us back to the original question:
how to represent the solution of the recurrence for *all points in simultaneous*?


### Propagating zonotopes

The set representation that is most effective to this problem are *zonotopes* since
they are closed under linear maps. Moreover, since the initial states is a
hyperrectangle (thus a zonotope), we can propagate the whole set exactly (and efficiently).

We now implement the solution by propagating zonotopes.

```@example discrete_propagation
# map X0 according to the rotation matrix
X0z = convert(Zonotope, X0)

arr = [linear_map(M(θi), X0z) for θi in range(0, 2pi, length=50)]

plot(arr, ratio=1., alpha=1.)
plot!(X0, lw=2.0, ls=:dash, alpha=1., c=:white)
```

```@example discrete_propagation
typeof(arr)
```

```@example discrete_propagation
X0
```

```@example discrete_propagation
convert(Zonotope, X0)
```

```@example discrete_propagation
genmat(arr[1])
```

(do `?Zonotope` for the docstring)

```@example discrete_propagation
genmat(arr[2])
```

```@example discrete_propagation
M(2pi/50) * genmat(arr[1])
```

The set effectively rotates clockwise around the origin:

```@example discrete_propagation
plot(Singleton(zeros(2)))
plot!(arr[1:10], ratio=1., alpha=1.)
plot!(X0, lw=2.0, ls=:dash, alpha=1., c=:white)
```

!!! note "Animations with Plots"
    The Julia package Plots.jl has a `@gif` functionality that can be used
    to created animated gifs. The following code animates the solution at 15 frames per second (fps).

    ```julia
    # map X0 according to the rotation matrix
    X0z = convert(Zonotope, X0)

    fig = plot()
    anim = @animate for θi in range(0, 2pi, length=50)
        plot!(fig, linear_map(M(θi), X0z), lw=2.0, alpha=1., lab="")
        xlims!(-3, 3)
        ylims!(-2, 2)
    end
    gif(anim, fps=15)
    ```

Now we will use *reach-sets* and associate the angle $\theta$ with the time field.

```@example discrete_propagation
# propagate sets passing the time point (angle)
Rsets = [ReachSet(M(θi) * X0, θi) for θi in range(0, 2pi, length=50)]
nothing # hide
```

We can pass an array of reach-sets to the plotting function:

```@example discrete_propagation
plot(Rsets, vars=(1, 2), xlab="x", ylab="y", ratio=1., c=:blue)
plot!(X0, c=:white, alpha=.6)
```

Since reach-sets have time information, we can also plot the sequence *in time*.

```@example discrete_propagation
plot(Rsets, vars=(0, 1), xlab="t", lab="x(t)", lw=2.0, lc=:blue, alpha=1.)
plot!(Rsets, vars=(0, 2), lab="y(t)", lw=2.0, lc=:orange, alpha=1.)
```

## What is a flowpipe?

A **flowpipe** is a collection of reach-sets which behaves like their (set) union. Flowpipes attain the right level of abstraction in order to represent solutions of set-based problems.

We can instantiate a flowpipe by passing an array of reach-sets.

```@example discrete_propagation
arr = [ReachSet(M(θi) * X0, θi) for θi in range(0, 2pi, length=50)]

F = Flowpipe(arr)
typeof(F)
```

We can plot flowpipes, and all the reach-sets are plotted with the same color.

```@example discrete_propagation
plot(F, vars=(1, 2), ratio=1.)
```

Flowpipes implement Julia's array inteface.

```@example discrete_propagation
length(F)
```

For instance, do `F[1:3:end]` to plot one every three elements:

```@example discrete_propagation
plot!(F[1:3:end], vars=(1, 2), c=:red)
```

Of course, it is also possible to use the wrapped array (do `array(F)`) directly. However, flowpipes can be used to filter reach-sets in time, among other operations.

```@example discrete_propagation
tspan.(F[1:5]) # time span of the first five reach-sets
```

```@example discrete_propagation
F(2pi) # find the reach-set at time t = 2pi
```

We can also pick, or filter, those reach-sets whose intersection is non-empty with
a given time interval.

```@example discrete_propagation
F(2.5 .. 3.0) # returns a view
```

We can also filter by a condition on the time span using `tstart` and `tend`:

```@example discrete_propagation
# get all reach-sets from time t = 4 onwards
aux = F(4 .. tend(F));

length(aux)
```

```@example discrete_propagation
plot(F, vars=(0, 1), xlab="t", ylab="x", lw=2.0, alpha=1.)

# get all those reach-set whose time span is greater than 40
plot!(aux, vars=(0, 1), lw=2.0, lc=:red, alpha=1.)
```

Finally, observe that set operations with flowpipes are also supported. The following
example intersects the flowpipe $F$ with a half-space.

```@example discrete_propagation
H  = HalfSpace([1, 1.], 1.) # x + y ≤ 1

Q = F ∩ H  # perform a lazy intersection

# plot the result
plot(H, alpha=.3, lab="H", c=:grey)
plot!(F ∩ H, vars=(1, 2), ratio=1., lab="F ∩ H")
xlims!(-2.0, 2.0); ylims!(-2, 2.)
```

### Using the solve interface for discrete problems

By default, the `solve` interface propagates sets in dense time (as explained in the next section).
However, it is possible to propagate sets in discrete time by specifying that the approximation model should not add any bloating to the initial states, `approx_model=NoBloating()` to the solver used.

**Exercise.** Solve the linear set based recurrence in the current section using the `solve` interface.
Hint: use the algorithm `alg=GLGM06(δ=0.05, approx_model=NoBloating())` and the state matrix as
specified in the following section.


## Set propagation in dense time

In this section we transition from the discrete time problem considered before
to a problem in which time is assumed to be a continuously varying quantity.
Thus, our next "toy" problem is the system of ordinary differential equations
that describe a simple, undamped, harmonic oscillator:

```math
\left\{ \begin{aligned}
x'(t) &= y(t) \\
y'(t) &= -x(t)
\end{aligned} \right. \qquad \textrm{subject to } (x(0), y(0)) \in X_0 = [0.8, 1.2] \times [0.8, 1.2] \subset \mathbb{R}^2
```
We are interested in solving this problem for all $t \in [0, T]$ with $T = 2\pi$.

The problem studied in the previous section is the "discrete analogue" of the purely continuous system of differential equations defined above.

The matrix $M(\omega t)$ describes the solution of a mathematical model called *simple harmonic oscillator* with natural frequency $\omega = 1$, and the **analytic solution** is

```math
\left\{ \begin{aligned}
x(t) &= x_0 \cos~ t + y_0 \sin~ t \\
y(t) &= -x_0 \sin t + y_0 \cos t
\end{aligned} \right.
```

Observe that $(x(t), y(t))^T$ is just the matrix $M(t)$ applied to the initial state $(x_0, y_0)^T$.


Define the invariant $G: \{(x, y) \in \mathbb{R}^2: x ≥ 1.3 \}$.


```@example dense_propagation
using ReachabilityAnalysis # hide
using ReachabilityAnalysis: center  # hide

import Plots: plot, plot!, xlims!, ylims! # hide

import Random # hide
Random.seed!(1117) # hide

# initial states
X0 = BallInf(ones(2), 0.2)

# rotation matrix
M(θ) = [cos(θ) sin(θ); -sin(θ) cos(θ)]

# analytic solution to visualize the solution at intermediate times
function analytic_sol(u0, tt)  # hide
    sol = [M(t) * u0 for t in tt]  # hide
    xcoords = getindex.(sol, 1)  # hide
    ycoords = getindex.(sol, 2)  # hide
    return xcoords, ycoords  # hide
end  # hide

samples = sample(X0, 200, include_vertices=true)
tt = range(0, 2pi/20, length=100)

F = [ReachSet(linear_map(M(ti), X0), ti) for ti in range(0, 2pi, step=2pi/20)]

fig = plot(F[1], vars=(1, 2), lab="X0", c=:blue)
plot!(fig, F[2], vars=(1, 2), lab="Φ * X0", c=:orange)
[plot!(fig, xcoords, ycoords, seriestype=:path, lab="", c=:magenta, ratio=1.) for (xcoords, ycoords) in analytic_sol.(samples, Ref(tt))]
fig
```

### Conservative time discretization



```@example dense_propagation
fig = plot()

plot!(fig, convexify(F[1:2]), vars=(1, 2), ls=:dash, lw=3.0, c=:green, lab="CH(X0, Φ*X0)")
plot!(fig, F[1], vars=(1, 2), lab="X0", c=:blue)
plot!(fig, F[2], vars=(1, 2), lab="Φ*X0", c=:orange)
[plot!(fig, xcoords, ycoords, seriestype=:path, lab="", c=:magenta, ratio=1.) for (xcoords, ycoords) in analytic_sol.(samples, Ref(tt))]

fig
```

```@example dense_propagation
xlims!(1.0, 1.5)
ylims!(0.8, 1.2)
```

```@raw html
<div style="background: ghostwhite;
            font-size: 16px;
            padding: 10px;
            border: 1px solid lightgray;
            margin: 10px;
            font-family: helvetica;">

<b>Proposition 1.</b> Conservative time discretization for linear homogeneous systems [[LGG09]](@ref).

Let $x' = Ax$ with initial states $x(0) \in X_0$. Then the convex set

$$
\Omega_0 := CH(X_0, \Phi X_0 \oplus B_r),
$$

where $\Phi = e^{A \delta}$ is the state transition matrix given a step size
$\delta > 0$ and $B_r$ is the infinity norm ball centered at the origin and
radius $r = (e^{\delta \Vert A \Vert } - 1 - \delta \Vert A \Vert)\Vert X_0\Vert$,
satisfies $R^e([0, \delta], X_0) \subseteq \Omega_0$, where
$R^e([0, \delta], X_0)$ are the exact reachable states between time $0$ and $\delta$.

Moreover, $\Omega_0$ converges to the true reachable states when the step size
$\delta \to 0$ (in Hausdorff distance).

</div>
```

Algorithms implementing conservative time discretization can be used from the
[`discretize(ivp::IVP, δ, alg::AbstractApproximationModel)`](@ref) function.
Set-based conservative discretization of a continuous-time initial value problem
into a discrete-time problem.
This function receives three inputs: the initial value problem (``ivp`) for a
linear ODE in canonical form, (e.g. the system returned by `normalize`);
the step-size (`δ`), and the algorithm (`alg`) used to compute the approximation model.


```@example
subtypes(ReachabilityAnalysis.AbstractApproximationModel)
```

The output of a discretization is a new initial value problem of a discrete system.
Different approximation algorithms and their respective options are described
in the docstring of each method, e.g. [`Forward`](@ref).

Initial-value problems considered in this function are of the form

```math
x' = Ax(t) + u(t),\\qquad x(0) ∈ \\mathcal{X}_0,\\qquad (1)
```
and where ``u(t) ∈ U(k)`` add where ``\\{U(k)\\}_k`` is a sequence of sets of
non-deterministic inputs and ``\\mathcal{X}_0`` is the set of initial
states. Recall that this initial-value problem is called homogeneous whenever `U`
is the empty set. Other problems, e.g. ``x' = Ax(t) + Bu(t)`` can be brought
to the canonical form with the function [`normalize`](@ref).

The initial value problem returned by this function consists of a set discretized
(also called *bloated*) initial states ``Ω₀``, together with the coefficient matrix
``Φ = e^{Aδ}`` and a transformed sequence of inputs if ``U`` is non-empty.

Two main variations of this algorithm are considered: dense time case and
discrete time case.

- In the dense time case, the transformation is such that the trajectories
of the given continuous system are included in the computed flowpipe of the
discretized system. More precisely, given a step size ``δ`` and the system (1)
conservative set-based discretization function computes a set, ``Ω₀``, that
guarantees to contain all the trajectories of (1) starting at any ``x(0) ∈ \\mathcal{X}_0``
and for any input function that satisfies ``u(t) ∈ U(1)``, for any ``t ∈ [0, δ]``.
If ``U`` is time-varying, this function also discretizes the inputs for ``k ≥ 0``.

- In the discrete time case, there is no bloating of the initial states and the
input is assumed to remain constant between sampled times. Use the algorithm
`NoBloating()` for this setting. If ``U`` is time-varying, this function also discretizes
the inputs for ``k ≥ 0``.

There are algorithms to obatin such transformations, called *approximation models*
in the technical literature. For references to the original papers, see the
docstring of each concrete subtype of `AbstractApproximationModel`.
