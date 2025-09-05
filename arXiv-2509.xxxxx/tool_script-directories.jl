# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

data_directory = joinpath(@__DIR__, "data")
external_data_directory = joinpath(data_directory, "ext")
output_data_directory = joinpath(data_directory, "out")
dump_data_directory = joinpath(output_data_directory, "dump")

external_code_directory = joinpath(@__DIR__, "external_code")
BlackHawk_directory = joinpath(external_code_directory, "BlackHawk_v2.3")
BlackHawk_hack_directory = joinpath(external_code_directory, "BlackHawk_v2.3-hack")
BlackHawk_results_directory = joinpath(BlackHawk_directory, "results")

plot_directory = joinpath(@__DIR__, "plots")
