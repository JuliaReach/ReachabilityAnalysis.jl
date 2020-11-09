import Literate
using Literate: script, markdown, notebook

# source directory for model files written using Literate.jl format
source_dir = joinpath(@__DIR__, "..", "examples")

# target directory for the tests files (bare .jl files)
target_dir_jl = joinpath(@__DIR__, "..", "test", "models")

# target directory for the markdown files (used by Documenter)
target_dir_md = joinpath(@__DIR__, "src", "models")
mkpath(target_dir_md)

# model files in sub-directories of source_dir
MODELS = ["Brusselator/Brusselator.jl",
          "Building/Building.jl",
          "DuffingOscillator/DuffingOscillator.jl",
          "EMBrake/EMBrake.jl",
          "Heat3D/Heat3D.jl",
          "ISS/ISS.jl",
          "LaubLoomis/LaubLoomis.jl",
          "Lorenz/Lorenz.jl",
          "LotkaVolterra/LotkaVolterra.jl",
          #"LotkaVolterraTangential/LotkaVolterraTangential.jl",
          "OpAmp/OpAmp.jl",
          "SquareWaveOscillator/SquareWaveOscillator.jl",
          "Platoon/Platoon.jl",
          "ProductionDestruction/ProductionDestruction.jl",
          "Quadrotor/Quadrotor.jl",
          "SEIR/SEIR.jl",
          "Spacecraft/Spacecraft.jl",
          "TransmissionLine/TransmissionLine.jl",
          "VanDerPol/VanDerPol.jl"]

for file in MODELS
    source_path = abspath(joinpath(source_dir, file))
    text = script(source_path, target_dir_jl, credit=false)
    code = strip(read(text, String))
    mdpost(str) = replace(str, "@__CODE__" => code)
    markdown(source_path, target_dir_md, postprocess=mdpost, credit=false)
    notebook(source_path, target_dir_md, execute=true, credit=false)
end
