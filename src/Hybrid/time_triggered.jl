# Hybrid systems with time-triggered transitions

function solve_hybrid(S, X0, reset_map;
                      Tsample=0.75, max_jumps=100, ζ=0.005,
                      t_offset=0.0, alg=GLGM06(δ=0.0003))

    # resuelvo primer tramo
    prob = IVP(S, X0)
    sol = solve(prob, T=2Tsample, GLGM06(δ=δ))

    # preallocate output vector of flowpipes
    # TODO: usar un tipo eg. HybridFlowpipe
    FT = typeof(flowpipe(sol))
    out = Vector{FT}()

    t0 = 0.0
    for k in 1:max_jumps

        # agrego segmento
        aux = Vector(sol(0 .. Tsample-ζ))
        # porque sol devuelve un view, despues
        # voy a hacer un flowpipe de view
        push!(out, shift(Flowpipe(aux), t0))

        t0 += tstart(aux[end]) # DOUBLE CHECK! see picture below:
        # tstart       Ts-ζ          tend
        # [-------------|-------------]

        Xend = sol(Tsample-ζ .. Tsample+ζ) |> Vector |> Flowpipe |> Convexify |> set

        sol = solve(IVP(S, reset_map(Xend)), T=2Tsample)

    end

    return out
end
