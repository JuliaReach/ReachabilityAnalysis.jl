## A warning note

Users of the library should have in mind that parameter tuning is an
essential ingredient to the successful application of reachability analysis.
This is in part due to the "cutting-edge research" aspect
of the methods available, i.e. the methods haven't yet been "battle-tested".
Moreover, finding good algorithm heuristics -- specially for hybrid systems -- is a hard
problem by itself.

On the other hand, our goal in designing and building the tools around `JuliaReach` has
been to make the default settings already work, or at least, work reasonably well across
a broad range of models: linear, hybrid, nonlinear, and so on. In the future, we hope that
users won't have to know the details of how the different set-based reachability methods
work under-the-hood in order to get new results. Eventually, users should know how
to adjust *some* parameters associated to the algorithms in order to improve the precision,
to formulate and solve verification problems, etc. inasmuch as effectively applying
standard numerical integration requires some basic knowledge about stepsize control,
interpolation, error calculation, etc.
