```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Flowpipes

## Abstract types

```@docs
AbstractFlowpipe
AbstractReachSet
```

## Interface functions

```@docs
set
setrep
tstart
tend
tspan
```

## Reachable set

```@docs
ReachSet
```

## Sparse reachable set

```@docs
SparseReachSet
```

## Type-stable flowpipe

```@docs
Flowpipe
```

## Hybrid flowpipe

```@docs
HybridFlowpipe
```

!!! note

  Write some conversion functions from e.g. Abaqus file format.

```julia
for face in 1:nfaces(cell)
    if onboundary(cell, face) && (cellid(cell), face) ∈ getfaceset(grid, "Neumann Boundary")
        reinit!(facevalues, cell, face)
        for q_point in 1:getnquadpoints(facevalues)
            dΓ = getdetJdV(facevalues, q_point)
            for i in 1:getnbasefunctions(facevalues)
                δu = shape_value(facevalues, q_point, i)
                fe[i] += δu * b * dΓ
            end
        end
    end
end
```

We start by looping over all the faces of the cell, next we have to check if
this particular face is located on the boundary, and then also check that the
face is located on our face-set called `"Neumann Boundary"`. If we have determined
that the current face is indeed on the boundary and in our faceset, then we
reinitialize `facevalues` for this face, using [`discretize`](@ref). When `reinit!`ing
`facevalues` we also need to give the face number in addition to the cell.
Next we simply loop over the quadrature points of the face, and then loop over
all the test functions and assemble the contribution to the force vector.

!!! note "Examples"
    The following commented examples makes use of Neumann boundary conditions:
    - TODO
