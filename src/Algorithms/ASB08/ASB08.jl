# pehaps call it CompoundFlowpipe / CompositeFlowpipe / PartitionedFlowpipe
# this a collection of flowpipes (each start time matches with the end time)
# but each element is itself a MixedFlowpipe, with potentially several Flowpipes
# This algorithm computes a sequence F₁, F₂, …, Fₖ of flowpipes such that
# with associated time spans (0, Δ), (2Δ, 3Δ), …, ((k-1)Δ, kΔ)

# running example: Van der Pol
n

N = Float64
ST = Zonotope{N, SVector{n, N}, SMatrix{n, n, N, n * n}}
FT = Flowpipe{N, ST, Vector{ST}}

HybridFlowpipe{}

# settings for the linear reachability algorithm
linearize and solve continuous system
δ
alg = GLGM06(δ=δ, static=true, dim=2, ngens=2, max_order=1)
tspan = 0 .. Δ

post(ivp, tspan, alg)
