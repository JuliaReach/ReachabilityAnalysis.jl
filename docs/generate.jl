import Literate
import ReachabilityBase.CurrentPath: @current_path

source_dir = joinpath(@__DIR__, "..", "examples")
# target directory for script files (*.jl) (used in the tests)
target_dir_jl = joinpath(@__DIR__, "..", "test", "models", "generated")
# target directory for markdown and notebook files (used by Documenter)
target_dir_md = joinpath(@__DIR__, "src", "generated_examples")
mkpath(target_dir_md)

MODELS = [
          #
          "Brusselator",
          "Building",
          "DuffingOscillator",
          "EMBrake",
          "ISS",
          "LaubLoomis",
          "Lorenz",
          "LotkaVolterra",
          "OpAmp",
          "SquareWaveOscillator",
          "Platoon",
          "ProductionDestruction",
          "Quadrotor",
          "SEIR",
          "Spacecraft",
          "TransmissionLine",
          "VanDerPol"
          #
          ]

macro current_path(prefix::String, filename::String)
    return joinpath(source_dir, prefix, filename)
end

for model in MODELS
    model_path = abspath(joinpath(source_dir, model))
    for file in readdir(model_path)
        if endswith(file, ".jl")
            input = abspath(joinpath(model_path, file))
            script = Literate.script(input, target_dir_jl; credit=false)
            code = strip(read(script, String))
            mdpost(str) = replace(str, "@__CODE__" => code)
            if get(ENV, "DOCUMENTATIONGENERATOR", "") == "true"
                Literate.markdown(input, target_dir_md; postprocess=mdpost, credit=false)
            else
                # for the local build, one needs to set `nbviewer_root_url`
                Literate.markdown(input, target_dir_md; postprocess=mdpost, credit=false,
                                  nbviewer_root_url="..")
            end

            # notebooks are deactivated to speed up the generation
            # Literate.notebook(input, target_dir_md; execute=true, credit=false)
            # if used, add the following to the top of the script files (where `MODELNAME` is the model name):
            #md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated_examples/MODELNAME.ipynb)
        elseif any(endswith.(file, [".jld2", ".png"]))
            # ignore *.jld2 and *.png files without warning
        else
            @warn "ignoring $file in $model_path"
        end
    end
end
