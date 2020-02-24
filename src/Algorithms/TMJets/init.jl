# out-of-place initialization
init(𝒫::TMJets, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝑆, copy(𝑂))

function options_TMJets()

    𝑂spec = Vector{OptionSpec}()

    # step size and termination criteria
    push!(𝑂spec, OptionSpec(:abs_tol, 1e-10, domain=Float64, info="absolute tolerance"))
    push!(𝑂spec, OptionSpec(:max_steps, 500, domain=Int, info="use this maximum number of steps in the validated integration"))

    # approximation options
    push!(𝑂spec, OptionSpec(:orderT, 10, domain=Int, info="order of the Taylor model in t"))
    push!(𝑂spec, OptionSpec(:orderQ, 2, domain=Int, info="order of the Taylor model for Jet transport variables"))

    # output options
    push!(𝑂spec, OptionSpec(:output_type, Hyperrectangle, info="output type of the Taylor model overapproximation"))

    return 𝑂spec
end

# in-place initialization
function init!(𝒫::TMJets, 𝑆::AbstractSystem, 𝑂::Options)

    # state dimension
    𝑂[:n] = statedim(𝑆)

    # adds default values for unspecified options
    𝑂init = validate_solver_options_and_add_default_values!(𝑂)

    return 𝑂init
end
