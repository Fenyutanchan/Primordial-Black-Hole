# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

import Pkg
import Pkg: PackageSpec

Pkg.activate(@__DIR__)

required_registered_packages = [
    "ArbNumerics",
    "CairoMakie",
    "CodecZlib",
    "DataInterpolations",
    "DelimitedFiles",
    "DifferentialEquations",
    "Format",
    "ForwardDiff",
    "IrrationalConstants",
    "JLD2",
    "LaTeXStrings",
    "PYTHIA8",
    "SpecialFunctions",
]

required_unregistered_packages = PackageSpec[
    PackageSpec(
        name="NaturalUnits",
        url="https://github.com/Fenyutanchan/NaturalUnits.jl.git"
    ),
]

try
    Pkg.add(required_unregistered_packages)
    Pkg.add(required_registered_packages)
catch
    rm(joinpath(@__DIR__, "Project.toml"); force=true, recursive=true)
    rm(joinpath(@__DIR__, "Manifest.toml"); force=true, recursive=true)
    Pkg.add(required_registered_packages)
    Pkg.add(required_unregistered_packages)
end
