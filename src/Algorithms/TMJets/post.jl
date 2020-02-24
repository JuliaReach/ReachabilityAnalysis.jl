using TaylorSeries: set_variables
using LazySets.Approximations: box_approximation

import TaylorModels
import IntervalArithmetic
const IA = IntervalArithmetic
const zeroI = IA.Interval(0.0) # [0.0, 0.0]
const oneI = IA.Interval(1.0) # [1.0, 1.0]
const symI = IA.Interval(-1.0, 1.0)
@inline zeroBox(m) = IntervalBox(zeroI, m)
@inline symBox(m) = IntervalBox(symI, m)

function _to_hyperrectangle(tTM, xTM, n)
    N = length(xTM)
    SET_TYPE = Hyperrectangle{Float64, SVector{n, Float64}, SVector{n, Float64}}
    Rsets = Vector{ReachSet{SET_TYPE}}(undef, N-1)
    @inbounds for i in 1:N-1
        Hi = convert(Hyperrectangle, xTM[i])
        t0 = tTM[i]
        t1 = tTM[i+1]
        Rsets[i] = ReachSet(Hi, t0, t1)
    end
    return Rsets
end

function _to_intervalbox(tTM, xTM, n)
    N = length(xTM)
    SET_TYPE = IA.IntervalBox{n, Float64}
    Rsets = Vector{ReachSet{SET_TYPE}}(undef, N-1)
    @inbounds for i in 1:N-1
        Bi = xTM[i]
        t0 = tTM[i]
        t1 = tTM[i+1]
        Rsets[i] = ReachSet(Bi, t0, t1)
    end
    return Rsets
end

function _to_zonotope(tTM, vTM, n)
    N = length(tTM)
    SET_TYPE = Zonotope{Float64}
    Rsets = Vector{ReachSet{SET_TYPE}}(undef, N-1)
    @inbounds for i in 1:N-1 # loop over the reach sets
        # pick the i-th Taylor model
        X = vTM[:, i]

        # pick the time domain of the given TM (same in all dimensions)
        Δt = TaylorModels.domain(X[1])

        # evaluate the Taylor model in time, the coefficents are now intervals
        X_Δt = evaluate(X, Δt)

        # builds the associated taylor model for each coordinate j = 1...n
        X̂ = [TaylorModelN(X_Δt[j], X[j].rem, zeroBox(n), symBox(n)) for j in 1:n]

        # floating point rigorous polynomial approximation
        fX̂ = TaylorModels.fp_rpa.(X̂)

        # LazySets can overapproximate a Taylor model with a Zonotope
        Zi = overapproximate(fX̂, Zonotope)
        t0 = tTM[i]
        t1 = tTM[i+1]

        Rsets[i] = ReachSet(Zi, t0, t1)
    end
    return Rsets
end

function post(𝒜::TMJets,
              𝑃::InitialValueProblem{<:Union{BBCS, CBBCS, CBBCCS}, <:LazySet},
              𝑂_global::Options)

    # ==================================
    # Initialization
    # ==================================

    𝑂 = merge(𝒜.options.defaults, 𝑂_global, 𝒜.options.specified)

    # system of ODEs
    f! = 𝑃.s.f
    n = statedim(𝑃)

    # initial time and time horizon
    t0 = 0.0
    T = 𝑂[:T]

    # maximum allowed number of steps
    max_steps = 𝑂[:max_steps]

    # unrap algorithm-specific options
    abs_tol, orderQ, orderT = 𝑂[:abs_tol], 𝑂[:orderQ], 𝑂[:orderT]

    # initial sets
    box_x0 = box_approximation(𝑃.x0)
    q0 = center(box_x0)
    δq0 = IntervalBox(low(box_x0)-q0, high(box_x0)-q0)

    # fix the working variables and maximum order in the global
    # parameters struct (_params_TaylorN_)
    set_variables("x", numvars=length(q0), order=2*orderQ)

    # define the property
    if 𝑂[:mode] == "check"
        property = 𝑂[:property]
    elseif 𝑂[:mode] == "reach"
        property = (t, x) -> true
    end

    # =====================
    # Flowpipe computation
    # =====================

    info("Reachable States Computation...")
    @timing begin
        tTM, xTM, vTM = validated_integ(f!, q0, δq0, t0, T, orderQ, orderT, abs_tol,
                                        maxsteps=max_steps, check_property=property)
    end

    # convert to hyperrectangle and wrap around the reach solution
    output_type = 𝑂[:output_type]
    if output_type == Hyperrectangle
        Rsets = _to_hyperrectangle(tTM, xTM, n)
    elseif output_type == IntervalBox
        Rsets = _to_intervalbox(tTM, xTM, n)
    elseif output_type == Zonotope
        Rsets = _to_zonotope(tTM, vTM, n)
    end

    Rsol = ReachSolution(Rsets, 𝑂)

    # ===========
    # Projection
    # ===========

    if 𝑂[:project_reachset]
        info("Projection...")
        Rsol = @timing project(Rsol)
    end

    return Rsol
end
