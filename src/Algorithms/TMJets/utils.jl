function _to_taylor_model(tTM, xTM, vTM, n)
    N = length(tTM)
    SET_TYPE = typeof(vTM[:, 1]) # this is Array{TaylorModel1{TaylorN{Float64},Float64},1}
    Rsets = Vector{ReachSet{SET_TYPE}}(undef, N-1)
    @inbounds for i in 1:N-1 # loop over the reach sets
        # pick the i-th Taylor model
        X = vTM[:, i]

        # pick the time domain of the given TM (same in all dimensions)
        Δt = TaylorModels.domain(X[1])

        Rsets[i] = ReachSet(X, IA.inf(Δt), IA.sup(Δt))
    end
    return Rsets
end

# Last edited: 27 feb 2020

# some constants
#const zeroI = IA.Interval(0.0) # [0.0, 0.0]
c#onst oneI = IA.Interval(1.0) # [1.0, 1.0]
#const symI = IA.Interval(-1.0, 1.0)
#@inline zeroBox(m) = IntervalBox(zeroI, m)
@inline unitBox(m) = IntervalBox(IA.Interval(0.0, 1.0), m)
#@inline symBox(m) = IntervalBox(symI, m)

# a taylor model reach set is a vector of taylor models (one for each state space variable)
# where each taylor model is the shape of a taylor model in one variable (time) whose coefficients
# are taylor models in several variables, ie. polynomials in the spatial variables
# moreover the spatial variables are normalized in -1 .. 1 (this is important!)
const TMReachSet = ReachSet{Vector{TaylorModel1{TaylorN{Float64},Float64}}}
const TMFlowpipe = Vector{TMReachSet}

# dimension of a taylor model reachset
function LazySets.dim(X::TMReachSet)
    return length(set(X)) # corresponds to the number of equations
end

# overapproximate taylor model reachset with a hyperrectangle
function LazySets.overapproximate(X::TMReachSet, ::Type{<:IntervalBox})
    # dimension of the reachset
    n = dim(X)

    # pick the time domain of the given TM (same in all dimensions)
    Δt = IA.Interval(X.t_start, X.t_end)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = evaluate(set(X), Δt)

    # evaluate the spatial variables in the symmetric box
    X̂ = IntervalBox([evaluate(X_Δt[i], symBox(n)) for i in 1:n]...)

    return X̂
end

function LazySets.overapproximate(X::TMReachSet, ::Type{<:IntervalBox}, nparts)
    # dimension of the reachset
    n = dim(X)

    # pick the time domain of the given TM (same in all dimensions)
    Δt = IA.Interval(X.t_start, X.t_end)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X_Δt = evaluate(set(X), Δt)

    # evaluate the spatial variables in the symmetric box
    partition = mince(symBox(n), nparts)
    X̂ = [IntervalBox([evaluate(X_Δt[i], Bi) for i in 1:n]...) for Bi in partition]

    return X̂
end
