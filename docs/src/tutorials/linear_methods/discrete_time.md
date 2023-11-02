```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
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
fig = plot(X0, c=:white)
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

import DisplayAs  # hide
fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
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

fig = plot(arr, ratio=1., alpha=1.)
plot!(X0, lw=2.0, ls=:dash, alpha=1., c=:white)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
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
fig = plot(Singleton(zeros(2)))
plot!(arr[1:10], ratio=1., alpha=1.)
plot!(X0, lw=2.0, ls=:dash, alpha=1., c=:white)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
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
fig = plot(Rsets, vars=(1, 2), xlab="x", ylab="y", ratio=1., c=:blue)
plot!(X0, c=:white, alpha=.6)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

Since reach-sets have time information, we can also plot the sequence *in time*.

```@example discrete_propagation
fig = plot(Rsets, vars=(0, 1), xlab="t", lab="x(t)", lw=2.0, lc=:blue, alpha=1.)
plot!(Rsets, vars=(0, 2), lab="y(t)", lw=2.0, lc=:orange, alpha=1.)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
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
fig = plot(F, vars=(1, 2), ratio=1.)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

Flowpipes implement Julia's array interface.

```@example discrete_propagation
length(F)
```

For instance, do `F[1:3:end]` to plot one every three elements:

```@example discrete_propagation
plot!(F[1:3:end], vars=(1, 2), c=:red)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
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
fig = plot(F, vars=(0, 1), xlab="t", ylab="x", lw=2.0, alpha=1.)

# get all those reach-set whose time span is greater than 40
plot!(aux, vars=(0, 1), lw=2.0, lc=:red, alpha=1.)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

Finally, observe that set operations with flowpipes are also supported. The following
example intersects the flowpipe $F$ with a half-space.

```@example discrete_propagation
H  = HalfSpace([1, 1.], 1.) # x + y ≤ 1

Q = F ∩ H  # perform a lazy intersection

# plot the result
fig = plot(H, alpha=.3, lab="H", c=:grey)
plot!(F ∩ H, vars=(1, 2), ratio=1., lab="F ∩ H")
xlims!(-2.0, 2.0); ylims!(-2, 2.)

fig = DisplayAs.Text(DisplayAs.PNG(fig))  # hide
```

### Using the solve interface for discrete problems

By default, the `solve` interface propagates sets in dense time (as explained in the next section).
However, it is possible to propagate sets in discrete time by specifying that the approximation model should not add any bloating to the initial states, `approx_model=NoBloating()` to the solver used.

**Exercise.** Solve the linear set based recurrence in the current section using the `solve` interface.
Hint: use the algorithm `alg=GLGM06(δ=0.05, approx_model=NoBloating())` and the state matrix as
specified in the following section.
