init(𝒜::ASB07, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒜, 𝑆, copy(𝑂))

function init!(::ASB07, 𝑆::AbstractSystem, 𝑂::Options)
    𝑂[:n] = statedim(𝑆)

    return validate_solver_options_and_add_default_values!(𝑂)
end
