# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function count_energy!(energy, energy_list, counter_list)
    @assert length(energy_list) == length(counter_list)
    ii_energy = findfirst(≥(energy), energy_list)

    if isnothing(ii_energy)
        counter_list[end] += 1
    else
        counter_list[ii_energy] += 1
    end

    return nothing
end
