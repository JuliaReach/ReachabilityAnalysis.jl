export ADB07d

struct ADB07d <: AbstractContinuousPost
    options::TwoLayerOptions

    function ADB07d(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_ASB07_decomposed())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
ADB07d(𝑂::Pair{Symbol,<:Any}...) =
    ADB07d(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
ADB07d() = ADB07d(Options())

include("options.jl")
include("init.jl")
include("post.jl")
include("reach.jl")
