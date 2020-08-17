import Literate
using Literate: script, markdown, notebook

source_dir = [joinpath(@__DIR__, "..", "examples")]
target_dir = joinpath(@__DIR__, "src", "models")
target_dir_test = joinpath(@__DIR__, "..", "test", "models")
mkpath(target_dir)

MODELS = ["VanDerPol/vanderpol.jl"]

for model in MODELS
    source_path = abspath(joinpath(source_dir, model))
    text = script(source_dir, target_dir_test, credit=false)
    code = strip(read(text, String))
    mdpost(str) = replace(str, "@__CODE__" => code)
    markdown(source_path, target_dir, postprocess=mdpost, credit=false)
    notebook(source_path, target_dir, execute=true, credit=false)
end
