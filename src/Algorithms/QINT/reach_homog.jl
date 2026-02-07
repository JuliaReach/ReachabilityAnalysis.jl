function _expansion_point(a, b, c, X0::Interval, Δ)
    c₀ = mid(X0.dat)
    f_c₀ = c₀ * (a * c₀ + b) + c
    x̃ = c₀ + Δ / 2 * f_c₀
    return x̃
end

function _linearize(a, b, c, x̃)
    μ = a * x̃
    α = 2 * μ + b
    β = -μ * x̃ + c
    return α, β
end

function reach_homog_QINT(; a, b, c, # right-hand side: f(x) = ax^2 + bx + c
                          X0,      # initial set
                          T,       # time span: [0, T]
                          Δ::N,    # step size for NL reach
                          δ::N,    # step size for the linear reach
                          θ,       # remainder expansion rate
                          maxiter) where {N}

    # total flowpipe
    RT = ReachSet{N,Interval{N}}
    VRT = Vector{RT}
    FT = Flowpipe{N,RT,VRT}
    Ftot = Vector{FT}()

    # initialization
    Δti = zeroI # defines the initial time
    waiting_list = WaitingList([Δti], [StateInLocation(X0, 1)])
    R̄err = Interval(-θ * Δ, θ * Δ)

    k = 1
    ksplit = 0
    while !isempty(waiting_list)
        (Δti, elem) = pop!(waiting_list)
        X0 = state(elem)

        # compute expansion point
        x̃ = _expansion_point(a, b, c, X0, Δ)

        # linearization
        α, β = _linearize(a, b, c, x̃)

        # solve linear reachability
        prob = @ivp(x' = α * x + β, x(0) ∈ X0)
        sol = post(INT(; δ=δ), prob, TimeIntervalC(0.0, Δ); Δt0=Δti)

        # compute admissible linearization error
        kθ = 1 / (exp(α * Δ) - 1) * α * θ * Δ
        L̄ = Interval(-kθ, kθ)

        # estimate lagrange remainder TODO: optimize
        R̂lin = ConvexHullArray([set(R) for R in sol]) ⊕ R̄err
        R̂lin_int = overapproximate(R̂lin, Interval)
        L = Interval(a * (R̂lin_int.dat - x̃)^2)

        # validate Lagrange remainder
        if L ⊆ L̄
            # if the inclusion check succeeded, add errors to each reach-set for final storage
            F = Flowpipe([ReachSet(Interval(set(R).dat + R̄err.dat), tspan(R)) for R in sol])
            push!(Ftot, F)

            # only add the final reach-set to the waiting list if we didn't hit the time horizon
            Rf = F[end]
            tf = tend(Rf)
            T <= tstart(Rf) && continue # TODO check
            elem = StateInLocation(set(Rf), 1)
            push!(waiting_list, tspan(Rf), elem)

        else
            ksplit += 1
            # if the inclusion didn't succeed, split initial condition
            # and append left and right interval segments to the waiting list
            X0l, X0r = split(X0, 2)
            iteml = StateInLocation(X0l, 1)
            itemr = StateInLocation(X0r, 1)

            push!(waiting_list, Δti, iteml)
            push!(waiting_list, Δti, itemr)
        end

        k += 1
        if k > maxiter
            error("maximum number of iterations maxiter=$maxiter reached, try increasing `maxiter`")
        end
    end # while

    return HybridFlowpipe(Ftot)
end
