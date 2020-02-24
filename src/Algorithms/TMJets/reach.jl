#=
Copyright (c) 2018-2019: David Sanders and Luis Benet.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

import IntervalArithmetic
using TaylorSeries
import TaylorIntegration
using TaylorModels: validated_step!, TaylorModel1, TaylorModelN

import IntervalArithmetic
const IA = IntervalArithmetic
const TI = TaylorIntegration

function validated_integ(f!, qq0::AbstractVector{T}, δq0::IntervalBox{N,T},
                         t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T,
                         params=nothing; # not currently used
                         maxsteps::Int=500, parse_eqs::Bool=true,
                         check_property::Function=(t, x)->true) where {N, T<:Real}

    # Set proper parameters for jet transport
    @assert N == get_numvars()
    dof = N

    # Some variables
    zI = zero(IA.Interval{T})
    q0 = IntervalBox(qq0)
    t = t0 + Taylor1(orderT)
    tI = t0 + Taylor1(orderT+1)
    δq_norm = IntervalBox(IA.Interval{T}(-1,1), N)
    q0box = q0 .+ δq_norm

    # Allocation of vectors
    # Output
    tv = Vector{T}(undef, maxsteps+1)
    xv = Vector{IntervalBox{N,T}}(undef, maxsteps+1)
    xTM1v = Matrix{TaylorModel1{TaylorN{T}, T}}(undef, dof, maxsteps+1)

    # Internals: jet transport integration
    x = Vector{Taylor1{TaylorN{T}}}(undef, dof)
    dx = Vector{Taylor1{TaylorN{T}}}(undef, dof)
    xaux = Vector{Taylor1{TaylorN{T}}}(undef, dof)
    xTMN = Vector{TaylorModelN{N, T, T}}(undef, dof)

    # Internals: Taylor1{Interval{T}} integration
    xI = Vector{Taylor1{IA.Interval{T}}}(undef, dof)
    dxI = Vector{Taylor1{IA.Interval{T}}}(undef, dof)
    xauxI = Vector{Taylor1{IA.Interval{T}}}(undef, dof)

    # Set initial conditions
    rem = Vector{IA.Interval{T}}(undef, dof)

    @inbounds for i in eachindex(x)
        qaux = normalize_taylor(qq0[i] + TaylorN(i, order=orderQ), δq0, true)
        x[i] = Taylor1(qaux, orderT)
        dx[i] = x[i]
        xTMN[i] = TaylorModelN(qaux, zI, q0, q0box)

        xI[i] = Taylor1(q0[i] + δq0[i], orderT+1)
        dxI[i] = xI[i]
        rem[i] = zI

        xTM1v[i, 1] = TaylorModel1(deepcopy(x[i]), zI, zI, zI)
    end

    # Output vectors
    @inbounds tv[1] = t0
    @inbounds xv[1] = evaluate(xTMN, δq_norm)

    # Determine if specialized jetcoeffs! method exists (built by @taylorize)
    parse_eqs = parse_eqs && (length(methods(TI.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TI.jetcoeffs!(Val(f!), t, x, dx, params)
        catch
            parse_eqs = false
        end
    end

    # Integration
    nsteps = 1
    while t0 < tmax

        # Validated step of the integration
        δt = validated_step!(f!, t, x, dx, xaux, tI, xI, dxI, xauxI,
            t0, tmax, xTMN, xv, rem, δq_norm, q0, q0box,
            nsteps, orderT, abstol, params, parse_eqs, check_property)

        # New initial conditions and time
        nsteps += 1
        t0 += δt
        @inbounds begin
            t[0] = t0
            tI[0] = t0
            tv[nsteps] = t0
            for i in eachindex(x)
                δtI = IA.Interval(0, δt)
                xTM1v[i, nsteps] = TaylorModel1(deepcopy(x[i]), rem[i], zI, δtI)
                aux = x[i](δt)
                x[i]  = Taylor1(aux, orderT)
                dx[i] = Taylor1(zero(aux), orderT)
                auxI = xTMN[i](δq_norm)
                xI[i] = Taylor1(auxI, orderT+1)
                dxI[i] = xI[i]
            end
        end

        if nsteps > maxsteps
            warn("Maximum number of integration steps $maxsteps reached; exiting. Try increasing `max_steps`.")
            break
        end

    end

    return view(tv, 1:nsteps), view(xv, 1:nsteps), view(xTM1v, :, 1:nsteps)
end
