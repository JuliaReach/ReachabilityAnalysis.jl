export ASB07

struct ASB07 <: AbstractContinuousPost
    options::TwoLayerOptions

    function ASB07(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_ASB07())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
ASB07(𝑂::Pair{Symbol,<:Any}...) = ASB07(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
ASB07() = ASB07(Options())

include("options.jl")
include("init.jl")
include("post.jl")
include("reach.jl")
