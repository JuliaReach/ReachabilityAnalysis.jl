# shared functionality among TMJets implementations

using TaylorModels: TaylorModelN
using TaylorModels: fp_rpa, remainder, initialize!

# =================================
# Defaut values for the parameters
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

# =================================================================
# Initialization funtions to prepare the input for validated_integ
# =================================================================

# fallback
function _initialize(X0::LazySet, orderQ, orderT)
    return _initialize(box_approximation(X0), orderQ, orderT)
end

# taylor model representations
_initialize(X0::TaylorModelReachSet, orderQ, orderT) = set(X0)
_initialize(X0::Vector{TaylorModel1{TaylorN{T}, T}}, orderQ, orderT) where {T} = X0

# hyperrectangular sets
_initialize(X0::AbstractHyperrectangle, orderQ, orderT) = convert(IntervalBox, box_approximation(X0))
_initialize(X0::IntervalBox, orderQ, orderT) = X0
_initialize(X0::IntervalArithmetic.Interval, orderQ, orderT) = IntervalBox(X0)

# zonotopic sets
function _initialize(X0::AbstractZonotope, orderQ, orderT)
    println("zonotopic set??")
    X = overapproximate(X0, TaylorModelReachSet, orderQ=orderQ, orderT=orderT)
    return set(X)
end

function _initialize(X0::CartesianProduct{N, <:AbstractZonotope, <:AbstractZonotope}, orderQ, orderT) where {N}
    println("USING THIS")
    X0z = convert(Zonotope, X0)
    return _initialize(X0z, orderQ, orderT)
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
    solver_kwargs = get(kwargs, :solver_kwargs, Dict(:maxsteps=>maxsteps))

    # call external solver
    tv, xv, xTM1v = solver_name(f!, X0box, t0, T, orderQ, orderT, abstol; parse_eqs, solver_kwargs...)

    # build flowpipe
    F = Vector{TaylorModelReachSet{N}}()
    sizehint!(F, maxsteps)
    for i in 2:length(tv)
        δt = TimeInterval(tv[i-1], tv[i])
        Ri = TaylorModelReachSet(xTM1v[:, i], δt + Δt0)
        push!(F, Ri)
    end
    ext = Dict{Symbol, Any}(:tv => tv, :xv => xv, :xTM1v => xTM1v)
    return Flowpipe(F, ext)
end
