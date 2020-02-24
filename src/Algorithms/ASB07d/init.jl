init(𝒜::ASB07_decomposed, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒜, 𝑆, copy(𝑂))

function init!(𝒜::ASB07_decomposed, 𝑆::AbstractSystem, 𝑂::Options)
    𝑂copy = copy(𝑂)
    𝑂copy[:n] = statedim(𝑆)

    𝑂validated = validate_solver_options_and_add_default_values!(𝑂copy)

    𝑂validated[:partition] = 𝒜.options[:partition]

    # compute blocks
    𝑂validated[:blocks] =
        compute_blocks(𝒜.options[:vars], 𝑂validated[:partition])

    return 𝑂validated
end
