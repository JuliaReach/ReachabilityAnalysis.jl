export solve

# ==============================================================================
# Algorithm defaults
# ==============================================================================

"""
    default_continuous_post(problem::InitialValueProblem)

Return the default continous post operator for the initial value problem of a
discrete or continuous system.

### Input

- `problem` -- an initial value problem represented by a mathematical system
               and a set of initial states

### Output

A continuous post operator with default options.
"""
function default_continuous_post(problem::InitialValueProblem{ST, XT}) where {ST<:AbstractContinuousSystem, XT}
    if islinear(problem)
        opC = GLGM06() # TODO consider: BFFPSV18()
    else
        opC = TMJets()
    end
    return opC
end

function default_discrete_post(problem::InitialValueProblem{ST, XT}) where {ST<:HybridSystem, XT}
    # TODO
end

function solve(problem::InitialValueProblem{ST, XT},
               opC::AbstractContinuousPost=default_continuous_post(problem);
               T::Float64,
               tspan::Tuple{Float64, Float64}) where {ST<:AbstractContinuousSystem, XT}

    # check arguments
    @assert statedim(problem) == dim(problem.x0) "the state-space dimension should match the " *
    "dimension of the initial states, but they are of size $(statedim(problem)) and $(dim(problem.x0)) respectively"
    # extract options
    # it is assumed  that either T or tspan are specified
    if !haskey(kwargs, :T)
        error("the time horizon `T` needs to be specified")
    end
    T = kwargs[:T]
    # consider tspan keyword argument

    # solve
    _solve_continuous(problem, opC; T=T)
end

# TODO hy hybrid system ?
function solve(problem::InitialValueProblem{ST, XT},
               opC::AbstractContinuousPost=default_continuous_post(problem),
               opD::AbstractContinuousPost=default_discrete_post(problem);
               T::Float64) where {ST<:HybridSystem, XT}

    # check arguments

    # solve
    _solve_hybrid(problem, opC, opD; T=T)
end

function _solve_continuous(problem, op; T...)

    # normalize system to canonical form if needed
    # TODO check normalization... probably only applies locally
    problem = IVP(normalize(problem.s), problem.x0)

    # run the continuous-post operator
    sol = post(problem, op; T=T)

    return sol
end

function _solve_hybrid(problem, post, kwargs)
    # TODO
end

#=

# ==============================================================================
# Main `solve` function for different system types
# ==============================================================================



## OLDER METHODS

"""
    solve(system::InitialValueProblem{<:HybridSystem},
          options::Options)

Interface to reachability algorithms for a hybrid system PWA dynamics.

### Input

- `system`  -- hybrid system
- `options` -- options for solving the problem
"""
function solve(system::InitialValueProblem{<:HybridSystem, <:LazySet},
               options::Options,
               opC::AbstractContinuousPost=BFFPSV18(),
               opD::AbstractDiscretePost=LazyDiscretePost())
    return solve!(distribute_initial_set(system), copy(options), opC, opD)
end

function solve!(system::InitialValueProblem{<:HybridSystem,
                                            <:Vector{<:Tuple{Int64,<:LazySet{N}}}},
               options_input::Options,
               opC::AbstractContinuousPost,
               opD::AbstractDiscretePost
              ) where N<:Real
    # update global variable
    global discrete_post_operator = opD

    HS = system.s
    init_sets = system.x0
    options = init!(opD, HS, options_input)
    time_horizon = options[:T]
    max_jumps = options[:max_jumps]
    property = options[:mode] == "check" ? options[:property] : nothing

    # waiting_list entries:
    # - (discrete) location
    # - (set of) continuous-time reach sets
    # - number of previous jumps
    waiting_list = Vector{Tuple{Int, AbstractReachSet{LazySet{N}}, Int}}()

    for (loc_id, x0) in init_sets
        loc = HS.modes[loc_id]
        source_invariant = loc.X

        if x0 ⊆ source_invariant
            loc_x0set = x0
        else
            loc_x0set = intersection(source_invariant, x0)
        end

        if !isempty(loc_x0set)
            aux = ReachSet{LazySet{N}}(loc_x0set, zero(N), zero(N))
            push!(waiting_list, (loc_id, aux, 0))
        end
    end

    # passed_list maps the (discrete) location to the (set of) continuous-time
    # reach sets
    passed_list = options[:fixpoint_check] == :none ?
        nothing :
        Vector{Vector{AbstractReachSet{LazySet{N}}}}(undef, nstates(HS))

    Rsets = Vector{AbstractReachSet{<:LazySet{N}}}()
    while (!isempty(waiting_list))
        loc_id, X0, jumps = pop!(waiting_list)
        loc = HS.modes[loc_id]
        source_invariant = loc.X

        # compute reach tube
        options_copy = copy(options)
        options_copy.dict[:T] = time_horizon - time_start(X0)
        options_copy.dict[:project_reachset] = false
        if property != nothing
            # TODO temporary hack, to be resolved in #447
            options_copy[:mode] = "reach"
        end

        if opC isa BFFPS19
            opC.options.specified[:HS] = HS
            opC.options.specified[:loc_id] = loc_id
            opC.options.specified[:opD] = opD
        end


        reach_tube = solve!(IVP(loc, set(X0)), options_copy, op=opC)

        # get the property for the current location
        property_loc = property isa Dict ?
                       get(property, loc_id, nothing) :
                       property

        # add the very first initial approximation
        if passed_list != nothing &&
                (!isassigned(passed_list, loc_id) || isempty(passed_list[loc_id]))
            reach_set = reach_tube.Xk[1]
            # TODO For lazy X0 the fixpoint check is likely to fail, so we
            # currently ignore that. In general, we want to add an
            # *underapproximation*, which we currently do not support.
            X = set(reach_set)
            if !(X isa CartesianProductArray) || !(array(X)[1] isa CH)
                Xoa = overapproximate(X)
                ti, tf = time_start(reach_set), time_end(reach_set)
                passed_list[loc_id] = [ReachSet{LazySet{N}}(Xoa, ti, tf)]
            end
        end

        # count_Rsets counts the number of new reach sets added to Rsets
        count_Rsets = tube⋂inv!(opD, reach_tube.Xk, loc.X, Rsets,
                                [time_start(X0), time_end(X0)])

        if property_loc != nothing
            if opD isa DecomposedDiscretePost
                temp_vars = opD.options[:temp_vars]
                n_lowdim = length(temp_vars)
                n = dim(set(X0))
                property_loc_lowdim = project(property_loc, temp_vars)
                for (i, reach_set) in enumerate(reach_tube.Xk)
                    X = set(reach_set)
                    if (dim(X) == n_lowdim && n_lowdim < n)
                        if !check(property_loc_lowdim, X)
                            return CheckSolution(false, i, options)
                        end
                    elseif !check(property_loc, X)
                        return CheckSolution(false, i, options)
                    end
                end
            else
                for (i, reach_set) in enumerate(Rsets[length(Rsets) - count_Rsets + 1 : end])
                    if !check(property_loc, set(reach_set))
                        return CheckSolution(false, i, options)
                    end
                end
            end
        end

        if jumps == max_jumps
            continue
        end
        post(opD, HS, waiting_list, passed_list, loc_id, Rsets, count_Rsets,
            jumps, options)

    end
    if options[:mode] == "check"
        return CheckSolution(true, -1, options)
    end

    # create vector with concrete set type (needed by ReachSolution)
    Rsets = [rs for rs in Rsets]

    # Projection
    if options[:project_reachset] || options[:projection_matrix] != nothing
        info("Projection...")
        RsetsProj = @timing project(Rsets, options)
    else
        RsetsProj = Rsets
    end
    return ReachSolution(RsetsProj, options)
end


import LazySets.use_precise_ρ

function postprocess(𝒫,
                     HS,
                     post_jump,
                     options,
                     waiting_list,
                     passed_list,
                     target_loc_id,
                     jumps)
    fixpoint_strategy = options[:fixpoint_check]

    if fixpoint_strategy == :eager
        # eager fixpoint checking
        post_jump_filtered =
            filter(x -> !isfixpoint(𝒫, x, passed_list, target_loc_id),
                   post_jump)
    else
        post_jump_filtered = post_jump
    end

    if (isempty(post_jump_filtered))
        # fixpoint found or transition can never be taken
        return
    end

    # apply clustering
    clustered = cluster(𝒫, post_jump_filtered, options)

    # push new sets after jump (unless a fixpoint is detected)
    for reach_set in clustered
        if fixpoint_strategy != :none
            if fixpoint_strategy == :lazy &&
                    isfixpoint(𝒫, reach_set, passed_list, target_loc_id)
                continue
            end
            push!(passed_list[target_loc_id], reach_set)
        end
        push!(waiting_list, (target_loc_id, reach_set, jumps))
    end
end

function cluster(𝒫::AbstractDiscretePost,
                 reach_sets::Vector{RSN},
                 options::Options) where {SN, RSN<:AbstractReachSet{SN}}
    clustering_strategy = options[:clustering]
    oa = 𝒫.options[:overapproximation]
    if clustering_strategy == :none
        # no clustering, keeping original set
        return reach_sets
    elseif clustering_strategy == :none_oa
        # no clustering but overapproximation
        return [RSN(overapproximate(set(reach_set), oa),
        time_start(reach_set), time_end(reach_set)) for reach_set in reach_sets]
    elseif clustering_strategy == :chull
        # cluster all sets in a convex hull and overapproximate that set
        chull = ConvexHullArray([set(reach_set) for reach_set in reach_sets])
        chull_oa = overapproximate(chull, oa)
        return [RSN(chull_oa, time_start(reach_sets[1]),
                    time_end(reach_sets[end]))]
    end
end

function isfixpoint(𝒫::AbstractDiscretePost,
                    reach_set::RSN,
                    passed_list,
                    loc_id
                   ) where {SN, RSN<:AbstractReachSet{SN}}
    @assert passed_list != nothing
    if isassigned(passed_list, loc_id)
        for other_reach_set in passed_list[loc_id]
            if set(reach_set) ⊆ set(other_reach_set)
                info("found a fixpoint in some reach tube")
                return true
            end
        end
        return false
    else
        passed_list[loc_id] = Vector{RSN}()
        return false
    end
end

# default: always apply line search
function use_precise_ρ(𝒫::AbstractDiscretePost,
                             cap::Intersection{N})::Bool where N<:Real
    return true
end

# --- default methods for handling assignments ---

# default implementation: use 'apply' from MathematicalSystems
function apply_assignment(𝒫::AbstractDiscretePost,
                          constrained_map::AbstractMap,
                          R⋂G::LazySet;
                          kwargs...)
    return apply(constrained_map, R⋂G)
end

# for reset maps: return a lazy ResetMap from LazySets
function apply_assignment(𝒫::AbstractDiscretePost,
                          constrained_map::ConstrainedResetMap,
                          R⋂G::LazySet;
                          kwargs...)
    return LazySets.ResetMap(R⋂G, constrained_map.dict)
end
=#
