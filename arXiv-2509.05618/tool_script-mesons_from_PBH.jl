# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

mesons_from_PBH = readdlm(joinpath(output_data_directory, "mesons_from_PBH.dat"); skipstart=1)
pion_from_PBH_interpolation = LinearInterpolation(log.(mesons_from_PBH[:, 2]), log.(mesons_from_PBH[:, 1]); extrapolation=ExtrapolationType.Linear, cache_parameters=true)
kaon_from_PBH_interpolation = LinearInterpolation(log.(mesons_from_PBH[:, 3]), log.(mesons_from_PBH[:, 1]); extrapolation=ExtrapolationType.Linear, cache_parameters=true)
kaon_zero_long_from_PBH_interpolation = LinearInterpolation(log.(mesons_from_PBH[:, 4]), log.(mesons_from_PBH[:, 1]); extrapolation=ExtrapolationType.Linear, cache_parameters=true)

dN_over_dt_pion(m_PBH) = (GeV ∘ exp ∘ pion_from_PBH_interpolation ∘ log ∘ EUval)(GeV, m_PBH)
dN_over_dt_kaon(m_PBH) = (GeV ∘ exp ∘ kaon_from_PBH_interpolation ∘ log ∘ EUval)(GeV, m_PBH)
dN_over_dt_kaon_zero_long(m_PBH) = (GeV ∘ exp ∘ kaon_zero_long_from_PBH_interpolation ∘ log ∘ EUval)(GeV, m_PBH)
