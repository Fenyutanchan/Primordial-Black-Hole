# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

Nₛ_data_from_PDG2024 = readdlm(joinpath(external_data_directory, "N_s_from_PDG2024.csv"), ',')
Nᵨ_OVER_Nₛ_data_from_PDG2024 = readdlm(joinpath(external_data_directory, "N_rho_OVER_N_s_from_PDG2024.csv"), ',')
N_data_around_mₑ = readdlm(joinpath(external_data_directory, "N_around_m_e.csv"), ',')

Nₛ_from_PDG2024_interpolation = LinearInterpolation(
    Nₛ_data_from_PDG2024[:, 2],
    log.(Nₛ_data_from_PDG2024[:, 1]);
    extrapolation=ExtrapolationType.Constant,
    cache_parameters=true
)
Nᵨ_OVER_Nₛ_from_PDG2024_interpolation = LinearInterpolation(
    Nᵨ_OVER_Nₛ_data_from_PDG2024[:, 2],
    log.(Nᵨ_OVER_Nₛ_data_from_PDG2024[:, 1]);
    extrapolation=ExtrapolationType.Constant,
    cache_parameters=true
)
Nᵨ_from_PDG2024_interpolation(log_T) = Nᵨ_OVER_Nₛ_from_PDG2024_interpolation(log_T) * Nₛ_from_PDG2024_interpolation(log_T)
Nᵨ_around_mₑ_interpolation = LinearInterpolation(
    N_data_around_mₑ[:, 2],
    log.(N_data_around_mₑ[:, 1]);
    extrapolation=ExtrapolationType.Constant,
    cache_parameters=true
)
Nₛ_around_mₑ_interpolation = LinearInterpolation(
    N_data_around_mₑ[:, 3],
    log.(N_data_around_mₑ[:, 1]);
    extrapolation=ExtrapolationType.Constant,
    cache_parameters=true
)

Nᵨ(EU::Type{<:EnergyUnit}, T::Real) = (Nᵨ ∘ EU)(T)
function Nᵨ(T::EnergyUnit)
    log_T_in_MeV = (log ∘ EUval)(MeV, T)
    T > MeV(10) && return Nᵨ_from_PDG2024_interpolation(log_T_in_MeV)
    return Nᵨ_around_mₑ_interpolation(log_T_in_MeV)
end
Nₛ(EU::Type{<:EnergyUnit}, T::Real) = (Nₛ ∘ EU)(T)
function Nₛ(T::EnergyUnit)
    log_T_in_MeV = (log ∘ EUval)(MeV, T)
    T > MeV(10) && return Nₛ_from_PDG2024_interpolation(log_T_in_MeV)
    return Nₛ_around_mₑ_interpolation(log_T_in_MeV)
end

function make_from_ρ_to_T(; Tₘᵢₙ_in_GeV=1e-5, Tₘₐₓ_in_GeV=1e3, number_of_points=1000)
    log10_T_in_GeV_list = log10.(geomspace(Tₘᵢₙ_in_GeV, Tₘₐₓ_in_GeV, number_of_points))
    log10_ρ_in_GeV_list = zero(log10_T_in_GeV_list)
    for (ii, log10_T_in_GeV) ∈ enumerate(log10_T_in_GeV_list)
        T = (GeV ∘ exp10)(log10_T_in_GeV)
        ρ = (π^2 / 30) * Nᵨ(T) * T^4
        log10_ρ_in_GeV_list[ii] = (log10 ∘ EUval)(GeV, ρ)
    end
    T_interpolation = LinearInterpolation(log10_T_in_GeV_list, log10_ρ_in_GeV_list; extrapolation=ExtrapolationType.Linear, cache_parameters=true)

    function from_ρ_to_T(ρ::EnergyUnit)
        @check_EU_dimension ρ 4

        ρ_val = EUval(GeV, ρ)
        ρ_val < 0 && return EU(NaN)
        ρ_val == 0 && return zero(EU)

        T_val = (exp10 ∘ T_interpolation ∘ log10)(ρ_val)

        return GeV(T_val)
    end

    return from_ρ_to_T
end
from_ρ_to_T = make_from_ρ_to_T()
