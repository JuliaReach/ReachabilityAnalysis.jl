# case without invariant: this is the same as in TaylorModels.jl but
# we wrap the each taylor model reach-set in the main loop
# TODO: do this a posteriri? then we can use `validated_integ!` directly
function _validated_integ!(F, f!, q0, δq0, t0, T, orderQ, orderT, abs_tol, max_steps, X::Universe, intersection_method)
    validated_integ!(F, f!, q0, δq0, t0, T, orderQ, orderT, abs_tol, max_steps)
end

# case with an invariant
function _validated_integ!(F, f!, q0, δq0, t0, T, orderQ, orderT, abs_tol, max_steps, X::LazySet, intersection_method)
    validated_integ!(F, f!, q0, δq0, t0, T, orderQ, orderT, abs_tol, max_steps, X, intersection_method)
end

# this function is the same as validated_integ! but with the addition that
# we consider an intersection with an invariant X in the main loop
# moreover we passs
function validated_integ!(F, f!, qq0::AbstractArray{T,1}, δq0::IntervalBox{N,T},
        t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T, max_steps::Int, X::LazySet,
        intersection_method::AbstractDisjointnessMethod, params=nothing;
        parse_eqs::Bool=true, check_property::Function=(t, x)->true) where {N, T<:Real}

    # Set proper parameters for jet transport
    @assert N == get_numvars()
    dof = N

    # Some variables
    zI = zero(IA.Interval{T})
    zbox = IntervalBox(zI, Val(N))
    symIbox = IntervalBox(IA.Interval{T}(-1, 1), Val(N))
    q0 = IntervalBox(qq0)
    t   = t0 + Taylor1(orderT)
    tI  = t0 + Taylor1(orderT+1)

    # Allocation of vectors
    # Output
    tv    = Array{T}(undef, max_steps+1)
    xv    = Array{IntervalBox{N,T}}(undef, max_steps+1)
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, max_steps+1)
    # Internals: jet transport integration
    x     = Array{Taylor1{TaylorN{T}}}(undef, dof)
    dx    = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xaux  = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xTMN  = Array{TaylorModelN{N,T,T}}(undef, dof)
    # Internals: Taylor1{IA.Interval{T}} integration
    xI    = Array{Taylor1{IA.Interval{T}}}(undef, dof)
    dxI   = Array{Taylor1{IA.Interval{T}}}(undef, dof)
    xauxI = Array{Taylor1{IA.Interval{T}}}(undef, dof)

    # Set initial conditions
    rem = Array{IA.Interval{T}}(undef, dof)
    @inbounds for i in eachindex(x)
        qaux = normalize_taylor(qq0[i] + TaylorN(i, order=orderQ), δq0, true)
        x[i] = Taylor1( qaux, orderT)
        dx[i] = x[i]
        xTMN[i] = TaylorModelN(qaux, zI, zbox, symIbox)
        #
        xI[i] = Taylor1( q0[i]+δq0[i], orderT+1 )
        dxI[i] = xI[i]
        rem[i] = zI
        #
        xTM1v[i, 1] = TaylorModel1(deepcopy(x[i]), zI, zI, zI)
    end
    sign_tstep = copysign(1, tmax-t0)

    # Output vectors
    @inbounds tv[1] = t0
    @inbounds xv[1] = TM.evaluate(xTMN, symIbox)

    # Determine if specialized jetcoeffs! method exists (built by @taylorize)
    parse_eqs = parse_eqs && (length(methods(TaylorIntegration.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(f!), t, x, dx, params)
        catch
            parse_eqs = false
        end
    end

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax

        # Validated step of the integration
        δt = validated_step!(f!, t, x, dx, xaux, tI, xI, dxI, xauxI,
            t0, tmax, sign_tstep, xTMN, xv, rem, zbox, symIbox,
            nsteps, orderT, abstol, params, parse_eqs, check_property)

        # New initial conditions and time
        nsteps += 1
        t0 += δt
        @inbounds t[0] = t0
        @inbounds tI[0] = t0
        @inbounds tv[nsteps] = t0
        @inbounds for i in eachindex(x)
            δtI = sign_tstep * IA.Interval{T}(0, sign_tstep*δt)
            xTM1v[i, nsteps] = TaylorModel1(deepcopy(x[i]), rem[i], zI, δtI)
            aux = x[i](δt)
            x[i]  = Taylor1( aux, orderT )
            dx[i] = Taylor1( zero(aux), orderT )
            auxI = xTMN[i](symIbox)
            xI[i] = Taylor1( auxI, orderT+1 )
            dxI[i] = xI[i]
        end

        # construct the taylor model reach-set
        Ri = TaylorModelReachSet(xTM1v[:, nsteps], TimeInterval(t0-δt, t0))

        # check intersction with invariant
        _is_intersection_empty(Ri, X) && break

        # update output flowpipe
        # note that F may have 1 less element than xTM1v and xv
        push!(F, Ri)

        if nsteps > max_steps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end

    end

    return F, view(tv,1:nsteps), view(xv,1:nsteps), view(xTM1v, :, 1:nsteps)
end
