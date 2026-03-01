# shared functionality among TMJets implementations

using TaylorModels: TaylorModelN, fp_rpa, remainder
using TaylorModels.ValidatedInteg: _DEF_MINABSTOL, validated_integ, validated_integ2

# =================================
# Default values for the parameters
# =================================

const DEFAULT_MAX_STEPS_TMJETS = 2000
const DEFAULT_ABS_TOL_TMJETS = 1e-10
const DEFAULT_ORDER_T_TMJETS = 8
const DEFAULT_ORDER_Q_TMJETS = 2

"""
    TMJets

The algorithm TMJets defaults to `TMJets21b`.
"""
const TMJets = TMJets21b

# =======================================================================
# Initialization functions to prepare the input for validated integration
# =======================================================================

# fallback
function _initialize(X0::LazySet, orderQ, orderT)
    return _initialize(box_approximation(X0), orderQ, orderT)
end

# Taylor model representations
function _initialize(X0::TaylorModelReachSet, orderQ, orderT)
    return set(X0)
end

_initialize(X0::Vector{TaylorModel1{TaylorN{T},T}}, orderQ, orderT) where {T} = X0

# hyperrectangular sets
function _initialize(X0::AbstractHyperrectangle, orderQ, orderT)
    IB = convert(IntervalBox, box_approximation(X0))
    return [IB[i] for i in 1:length(IB)]
end
_initialize(X0::IntervalBox, orderQ, orderT) = [X0[i] for i in 1:length(X0)]
_initialize(X0::IA.Interval, orderQ, orderT) = [X0]

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
    IB = convert(IntervalBox, convert(Hyperrectangle, X0))
    return [IB[i] for i in 1:length(IB)]
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
