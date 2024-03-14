# shared functionality among TMJets implementations

using TaylorModels: TaylorModelN
using TaylorModels: fp_rpa, remainder, initialize!

# =================================
# Default values for the parameters
# =================================

const DEFAULT_MAX_STEPS_TMJETS = 2000
const DEFAULT_ABS_TOL_TMJETS = 1e-10
const DEFAULT_ORDER_T_TMJETS = 8
const DEFAULT_ORDER_Q_TMJETS = 2

"""
    TMJets

The algorithm TMJets defaults to `TMJets2`.
"""
const TMJets = TMJets2

# =======================================================================
# Initialization functions to prepare the input for validated integration
# =======================================================================

# fallback
function _initialize(X0::LazySet, orderQ, orderT)
    return _initialize(box_approximation(X0), orderQ, orderT)
end

# taylor model representations
function _initialize(X0::TaylorModelReachSet, orderQ, orderT)
    return set(X0)
end

_initialize(X0::Vector{TaylorModel1{TaylorN{T},T}}, orderQ, orderT) where {T} = X0

# hyperrectangular sets
function _initialize(X0::AbstractHyperrectangle, orderQ, orderT)
    return convert(IntervalBox, box_approximation(X0))
end
_initialize(X0::IntervalBox, orderQ, orderT) = X0
_initialize(X0::IA.Interval, orderQ, orderT) = IntervalBox(X0)

# zonotopic sets
function _initialize(X0::AbstractZonotope, orderQ, orderT)
    X0z = convert(Zonotope, X0)
    if order(X0z) == 1
        X = overapproximate(X0z, TaylorModelReachSet; orderQ=orderQ, orderT=orderT)
    else
        X = overapproximate(X0z, TaylorModelReachSet; orderQ=orderQ, orderT=orderT,
                            box_reduction=true)
    end
    return set(X)
end

function _initialize(X0::CartesianProduct{N,<:AbstractHyperrectangle,<:AbstractHyperrectangle},
                     orderQ, orderT) where {N}
    return convert(IntervalBox, convert(Hyperrectangle, X0))
end

function _initialize(X0::CartesianProduct{N,<:Zonotope,<:Interval}, orderQ, orderT) where {N}
    Z = X0.X
    G = Z.generators
    ord = order(Z)
    n = dim(Z)

    if ord == 1
        # "exact"
        X = _overapproximate_structured(X0, TaylorModelReachSet; orderQ=orderQ, orderT=orderT)

    elseif (size(G) == (n, 2n + 1)) && isdiag(view(G, :, (n + 2):(2n + 1)))
        X = _overapproximate_structured_full(X0, TaylorModelReachSet; orderQ=orderQ, orderT=orderT)

    elseif (size(G) == (n, 2(n + 1))) && isdiag(view(G, :, (n + 2):(2n + 1))) &&
           iszero(view(G, :, 2n + 2))
        aux = Zonotope(Z.center, view(G, :, 1:(2n + 1)))
        X0z = convert(Zonotope, aux × X0.Y)
        X = _overapproximate_structured(X0z, TaylorModelReachSet; orderQ=orderQ, orderT=orderT)

    else # otherwise, resort to a box overapproximation
        X0z = convert(Zonotope, X0)
        X = overapproximate(X0z, TaylorModelReachSet; orderQ=orderQ, orderT=orderT,
                            box_reduction=true)
    end

    return set(X)
end

# =================================
# External solver function
# =================================

function _solve_external(f!, X0, t0, T, orderQ, orderT, abstol, maxsteps, Δt0; kwargs...)
    N = eltype(X0)

    # box overapproximation of the initial states
    X0box = convert(IntervalBox, box_approximation(X0))

    # extract solver name and options
    parse_eqs = get(kwargs, :parse_eqs, false)
    solver_name = get(kwargs, :solver_name, TM.validated_integ)
    solver_kwargs = get(kwargs, :solver_kwargs, Dict(:maxsteps => maxsteps))

    # call external solver
    tv, xv, xTM1v = solver_name(f!, X0box, t0, T, orderQ, orderT, abstol; parse_eqs,
                                solver_kwargs...)

    # build flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, maxsteps)
    @inbounds for i in eachindex(tv)
        δt = TimeInterval(tv[i] + domain(xTM1v[1, i]))
        Ri = TaylorModelReachSet(xTM1v[:, i], δt + Δt0)
        push!(F, Ri)
    end
    ext = Dict{Symbol,Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v)
    return Flowpipe(F, ext)
end
