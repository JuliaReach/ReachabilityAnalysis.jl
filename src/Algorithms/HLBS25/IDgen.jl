mutable struct IDGenerator
    max_id::Int
end

function synchronize!(idg::IDGenerator, xs...)
    all_ids = vcat(getfield.(xs, :idx)...)
    idg.max_id = maximum(all_ids)
    return idg
end

@inline function (idg::IDGenerator)(l::Int)
    @assert l ≥ 0 "number of IDs to generate must be ≥ 0"

    # prevent integer overflow
    if gen.max_id > typemax(Int) - l
        throw(OverflowError("ID space exhausted"))
    end

    start = gen.max_id + 1
    ids = Vector{Int}(undef, l)

    @inbounds @simd for i in 1:l
        ids[i] = start + i - 1
    end

    gen.max_id += l
    return ids
end

function fresh!(idg::IDGenerator,
                P::SparsePolynomialZonotope,
                As::MatrixZonotope...)

    idₚ = indexvector(P)
    # collect all IDs appearing in As
    idmz = Set{Int}(el for A in As for el in indexvector(A))

    # compute which IDs are shared
    shared = intersect(idmz, idₚ)

    # count how many need replacement
    n = 0
    @inbounds for id in idₚ
        n += !(id in shared)
    end

    # generate new unique IDs
    new_ids = idg(n)

    # replace
    j = 1
    @inbounds for i in eachindex(idₚ)
        if !(idₚ[i] in shared)
            idₚ[i] = new_ids[j]
            j += 1
        end
    end

    return P
end


function fresh!(idg::IDGenerator, P::SparsePolynomialZonotope)
    idₚ = indexvector(P)

    new_ids = idg(length(idₚ))

    @inbounds @simd for i in eachindex(idₚ)
        idₚ[i] = new_ids[i]
    end

    return P
end
