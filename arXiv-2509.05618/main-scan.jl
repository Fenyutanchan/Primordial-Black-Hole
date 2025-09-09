# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

# Preliminaries

using Distributed

# workers() == [1] || (rmprocs ∘ workers)()

# addprocs(24, exeflags="--project=$(@__DIR__)")
@info """
Added $(nworkers()) workers for $(basename(@__FILE__)).
Please run the following command to start the workers if the number of workers is not your desired one:
`julia --project --procs=<number of workers> $(basename(@__FILE__))`
"""

@everywhere using DelimitedFiles,
    Statistics

@everywhere using ArbNumerics,
    DataInterpolations, 
    DifferentialEquations,
    Format,
    ForwardDiff,
    JLD2,
    SpecialFunctions

@everywhere using NaturalUnits

@everywhere include("tool_script-ArbNumerics_compat.jl")
@everywhere include("tool_script-directories.jl")
@everywhere include("tool_script-geomspace.jl")
@everywhere include("tool_script-parameters.jl")

@everywhere include("tool_script-N_relativistic_species.jl")

@everywhere include("tool_script-mesons_from_PBH.jl")

@everywhere γ_PBH = inv(3 * sqrt(3))
@everywhere N_BH(T_BH) = (27 / 4) * Nᵨ(T_BH) / (30720 * π)
@everywhere T_BH(m_BH) = NU.M_Pl^2 / m_BH
@everywhere t_PBH_formation(m_PBH) = m_PBH / (γ_PBH * NU.m_Pl^2)

function evolution_system_SBBN!(du, u, p, η)
    a = exp(η)
    log_Xₙ, log_ρ_SM_in_EU = u

    Xₙ = exp(log_Xₙ)
    ρ_SM = EU(exp(log_ρ_SM_in_EU), 4)

    T_SM = from_ρ_to_T(ρ_SM)

    x = Q / T_SM
    Γₙ_SM = (255 / τ_n) * (12 + 6 * x + x^2) / x^5

    ### Xₙ evolution ##################################
    d_Xₙ_d_t_SM = -Γₙ_SM * (Xₙ - (1 - Xₙ) * exp(-x))
    d_Xₙ_d_t_decay = -Xₙ / τ_n

    d_Xₙ_d_t = d_Xₙ_d_t_SM + d_Xₙ_d_t_decay
    ##################################################

    ### energy density background and Hubble parameter
    H = sqrt(ρ_SM / (3 * NU.M_Pl^2))

    d_ln_N_ρ_d_ln_T = EUval(GeV, T_SM) * ForwardDiff.derivative(T_in_GeV -> (log ∘ Nᵨ)(GeV, T_in_GeV), EUval(GeV, T_SM))
    d_ln_N_s_d_ln_T = EUval(GeV, T_SM) * ForwardDiff.derivative(T_in_GeV -> (log ∘ Nₛ)(GeV, T_in_GeV), EUval(GeV, T_SM))

    d_log_ρ_SM_d_η = -4 + d_ln_N_ρ_d_ln_T / (3 + d_ln_N_s_d_ln_T)
    ##################################################

    du[1] = d_Xₙ_d_t / (H * Xₙ)
    du[2] = d_log_ρ_SM_d_η
    return du
end

evolution_problem_SBBN = ODEProblem(evolution_system_SBBN!,
    map(log ∘ (v -> EUval(EU, v)),
        [inv(1 + exp(Q / GeV(10))), (π^2 / 30) * Nᵨ(GeV(10)) * GeV(10)^4]
    ),
    (0, log(1e7))
)
evolution_solution_SBBN = solve(evolution_problem_SBBN, AutoTsit5(Rosenbrock23()))

a_list_SBBN = geomspace(1, 1e7, 1000)
results_list_SBBN = (evolution_solution_SBBN ∘ log).(a_list_SBBN)
ρ_SM_list_SBBN = map((x -> EU(x, 4)) ∘ exp ∘ last, results_list_SBBN)
T_SM_list_SBBN = from_ρ_to_T.(ρ_SM_list_SBBN)
Xₙ_list_SBBN = (exp ∘ first).(results_list_SBBN)

Yₚ_SBBN = begin
    tmp_perm = sortperm(Q ./ T_SM_list_SBBN)
    Xₙ_T_interpolation_SBBN = LinearInterpolation(Xₙ_list_SBBN[tmp_perm], Q ./ T_SM_list_SBBN[tmp_perm]; extrapolation=ExtrapolationType.Linear, cache_parameters=true)
    Xₙ_nuc_SBBN = Xₙ_T_interpolation_SBBN(Q / T_nuc)
    2 * Xₙ_nuc_SBBN / (1 - Xₙ_nuc_SBBN)
end

@info "SBBN scenario done with Yₚ = $(Yₚ_SBBN)."

@everywhere function evolution_system!(du, u, p, η)
    a = exp(η)
    log_t_in_EU, log_ρ_SM_in_EU, log_m_PBH_in_EU, log_Xₙ,
        log_n_pion_minus_in_EU, log_n_pion_plus_in_EU,
        log_n_kaon_minus_in_EU, log_n_kaon_zero_long_in_EU = u
    a_ini, n_PBH_ini, log_m_PBH_evaporated_in_EU = p

    PBH_evaporated_flag = log_m_PBH_in_EU ≤ log_m_PBH_evaporated_in_EU

    t = EU(exp(log_t_in_EU), -1)
    ρ_SM = EU(exp(log_ρ_SM_in_EU), 4)
    m_PBH = EU(exp(log_m_PBH_in_EU), 1)
    Xₙ = exp(log_Xₙ)
    n_pion_minus = EU(exp(log_n_pion_minus_in_EU < -100 ? -100 : log_n_pion_minus_in_EU), 3)
    n_pion_plus = EU(exp(log_n_pion_plus_in_EU < -100 ? -100 : log_n_pion_plus_in_EU), 3)
    n_kaon_minus = EU(exp(log_n_kaon_minus_in_EU < -100 ? -100 : log_n_kaon_minus_in_EU), 3)
    n_kaon_zero_long = EU(exp(log_n_kaon_zero_long_in_EU < -100 ? -100 : log_n_kaon_zero_long_in_EU), 3)

    n_PBH = PBH_evaporated_flag ? zero(n_PBH_ini) : n_PBH_ini * (a_ini / a)^3
    T_PBH = T_BH(m_PBH)
    T_SM = from_ρ_to_T(ρ_SM)

    n_γ = 2 * ζ₃ * T_SM^3 / π^2
    n_B = η_B * n_γ

    x = Q / T_SM
    ρ_PBH = m_PBH * n_PBH

    N_PBH = N_BH(T_PBH)
    Γₙ_SM = (255 / τ_n) * (12 + 6 * x + x^2) / x^5
    ρ_total = ρ_SM + ρ_PBH

    H = sqrt(ρ_total / (3 * NU.M_Pl^2))
    d_log_t_d_η = inv(H * t)

    d_ln_N_ρ_d_ln_T = EUval(GeV, T_SM) * ForwardDiff.derivative(T_in_GeV -> (log ∘ Nᵨ)(GeV, T_in_GeV), EUval(GeV, T_SM))
    d_ln_N_s_d_ln_T = EUval(GeV, T_SM) * ForwardDiff.derivative(T_in_GeV -> (log ∘ Nₛ)(GeV, T_in_GeV), EUval(GeV, T_SM))

    d_m_PBH_d_t = -N_PBH * NU.m_Pl^4 / m_PBH^2
    PBH_evaporated_flag && (d_m_PBH_d_t = zero(d_m_PBH_d_t))

    d_log_m_PBH_d_η = d_m_PBH_d_t / (H * m_PBH)
    # d_log_ρ_SM_d_η = -4 - d_m_PBH_d_t * n_PBH / (H * ρ_SM)
    d_log_ρ_SM_d_η = -4 + d_ln_N_ρ_d_ln_T / (3 + d_ln_N_s_d_ln_T) - d_m_PBH_d_t * n_PBH / (H * ρ_SM)

    ### Xₙ evolution ##################################
    d_Xₙ_d_t_SM = -Γₙ_SM * (Xₙ - (1 - Xₙ) * exp(-x))
    d_Xₙ_d_t_decay = -Xₙ / τ_n
    d_Xₙ_d_t_pion = (1 - Xₙ) * n_pion_minus * σv_pion_minus_pton(T_SM) - Xₙ * n_pion_plus * σv_pion_plus_ntop(T_SM)
    d_Xₙ_d_t_kaon_minus = (1 - Xₙ) * n_kaon_minus * σv_kaon_minus_pton(T_SM) - Xₙ * n_kaon_minus * σv_kaon_minus_ntop(T_SM)
    d_Xₙ_d_t_kaon_zero_long = (1 - Xₙ) * n_kaon_zero_long * σv_kaon_zero_long_pton(T_SM) - Xₙ * n_kaon_zero_long * σv_kaon_zero_long_ntop(T_SM)

    d_Xₙ_d_t = d_Xₙ_d_t_SM + d_Xₙ_d_t_decay +
                d_Xₙ_d_t_pion + d_Xₙ_d_t_kaon_minus +
                d_Xₙ_d_t_kaon_zero_long

    d_Xₙ_d_η = d_Xₙ_d_t / (H * Xₙ)
    ##################################################

    ### number density of mesons evolution ###########
    d_n_pion_d_t_from_PBH = dN_over_dt_pion(m_PBH) * n_PBH
    d_n_kaon_d_t_from_PBH = dN_over_dt_kaon(m_PBH) * n_PBH
    d_n_kaon_zero_long_d_t_from_PBH = dN_over_dt_kaon_zero_long(m_PBH) * n_PBH

    d_n_pion_minus_d_t = d_n_pion_d_t_from_PBH - Γ_pion_minus * n_pion_minus -
                            σv_pion_minus_pton(T_SM) * (1 - Xₙ) * n_B * n_pion_minus
    d_n_pion_plus_d_t = d_n_pion_d_t_from_PBH - Γ_pion_plus * n_pion_plus -
                            σv_pion_plus_ntop(T_SM) * Xₙ * n_B * n_pion_plus
    d_n_kaon_minus_d_t = d_n_kaon_d_t_from_PBH - Γ_kaon_minus * n_kaon_minus -
                            σv_kaon_minus_pton(T_SM) * (1 - Xₙ) * n_B * n_kaon_minus -
                            σv_kaon_minus_ntop(T_SM) * Xₙ * n_B * n_kaon_minus
    d_n_kaon_zero_long_d_t = d_n_kaon_zero_long_d_t_from_PBH - Γ_kaon_zero_long * n_kaon_zero_long -
                            σv_kaon_zero_long_pton(T_SM) * (1 - Xₙ) * n_B * n_kaon_zero_long -
                            σv_kaon_zero_long_ntop(T_SM) * Xₙ * n_B * n_kaon_zero_long
    ##################################################

    du[1] = d_log_t_d_η
    du[2] = d_log_ρ_SM_d_η
    du[3] = d_log_m_PBH_d_η
    du[4] = d_Xₙ_d_η
    du[5] = iszero(n_pion_minus) ? 0 : d_n_pion_minus_d_t / (H * n_pion_minus)
    du[6] = iszero(n_pion_plus) ? 0 : d_n_pion_plus_d_t / (H * n_pion_plus)
    du[7] = iszero(n_kaon_minus) ? 0 : d_n_kaon_minus_d_t / (H * n_kaon_minus)
    du[8] = iszero(n_kaon_zero_long) ? 0 : d_n_kaon_zero_long_d_t / (H * n_kaon_zero_long)

    return du
end

@everywhere function magic(m_PBH, β_PBH;
    BITS::Int=BITS, log_m_PBH_evaporated_ratio::Number=20, num_points::Int=10000,
    dump_flag::Bool=false, dump_num_points::Int=100,
    mesons_output_flag::Bool=false
)
    dump_flag && @assert dump_num_points <= num_points
    dump_indices = dump_flag ? (Int ∘ round).(1:((num_points - 1) / (dump_num_points - 1)):num_points) : nothing

    a_f = ArbFloat(1; bits=BITS)
    t_f = t_PBH_formation(m_PBH)
    H_f = inv(2 * t_f)
    ρ_SM_f = 3 * NU.M_Pl^2 * H_f^2
    T_SM_f = from_ρ_to_T(ρ_SM_f)
    ρ_PBH_f = β_PBH * ρ_SM_f
    n_PBH_f = ρ_PBH_f / m_PBH
    Xₙ_f = inv(1 + exp(Q / T_SM_f))

    a_end = 10 * cbrt(Nᵨ(T_SM_f) * a_f^3 * T_SM_f^3 / (Nᵨ(T_nuc) * T_nuc^3))

    n_meson_f = 1e-60 * 2 * ζ₃ * T_SM_f^3 / π^2

    # log_m_PBH_evaporated_in_EU = (log ∘ EUval)(EU, m_PBH) - log_m_PBH_evaporated_ratio
    log_m_PBH_evaporated_in_EU = max((log ∘ EUval)(EU, m_PBH) - log_m_PBH_evaporated_ratio, (log ∘ EUval)(EU, NU.m_Pl))

    evolution_problem = ODEProblem(evolution_system!,
        map(log ∘ (v -> EUval(EU, v)),
            [t_f, ρ_SM_f, m_PBH, Xₙ_f, n_meson_f, n_meson_f, n_meson_f, n_meson_f]
        ),
        (log(a_f), log(a_end)),
        (a_f, n_PBH_f, log_m_PBH_evaporated_in_EU)
    )
    evolution_solution = solve(evolution_problem, AutoTsit5(Rosenbrock23()))

    a_list = geomspace(a_f, a_end, num_points)
    results_list = (evolution_solution ∘ log).(a_list)

    t_list = map((v -> EU(v, -1)) ∘ exp ∘ (results -> results[1]), results_list)
    ρ_SM_list = map((v -> EU(v, 4)) ∘ exp ∘ (results -> results[2]), results_list)
    m_PBH_list = map((v -> EU(v, 1)) ∘ exp ∘ (results -> results[3]), results_list)
    n_PBH_list = n_PBH_f .* (a_f ./ a_list).^3
    Xₙ_list = map(exp ∘ (results -> results[4]), results_list)
    n_pion_minus_list = map((v -> EU(v, 3)) ∘ exp ∘ (results -> results[5]), results_list)
    n_pion_plus_list = map((v -> EU(v, 3)) ∘ exp ∘ (results -> results[6]), results_list)
    n_kaon_minus_list = map((v -> EU(v, 3)) ∘ exp ∘ (results -> results[7]), results_list)
    n_kaon_zero_long_list = map((v -> EU(v, 3)) ∘ exp ∘ (results -> results[8]), results_list)

    first_evaporated_index = findfirst(≤(EU(exp(log_m_PBH_evaporated_in_EU), 1)), m_PBH_list)
    isnothing(first_evaporated_index) || for ii ∈ first_evaporated_index:length(n_PBH_list)
        n_PBH_list[ii] = zero(n_PBH_list[ii])
    end

    # ρ_PBH_list = m_PBH_list .* n_PBH_list
    T_SM_list = from_ρ_to_T.(ρ_SM_list)

    # ρ_γ_list = (π^2 / 30) * 2 * T_SM_list.^4

    # ρ_total_list = ρ_SM_list + ρ_PBH_list

    Yₚ = begin
        tmp_perm = sortperm(Q ./ T_SM_list)
        Xₙ_T_interpolation_SBBN = LinearInterpolation(Xₙ_list[tmp_perm], Q ./ T_SM_list[tmp_perm]; extrapolation=ExtrapolationType.Linear, cache_parameters=true)
        Xₙ_nuc = Xₙ_T_interpolation_SBBN(Q / T_nuc)
        2 * Xₙ_nuc / (1 - Xₙ_nuc)
    end

    dump_flag && jldopen(joinpath(dump_data_directory, "dump_$(cfmt("%.2e", m_PBH / NU.g))gram_beta$(cfmt("%.2e", β_PBH)).jld2"), "w+"; compress=true) do io
        write(io, "a_list", Float64.(a_list[dump_indices]))
        write(io, "t_list [sec]", Float64.(t_list[dump_indices] ./ NU.s))
        write(io, "ρ_SM_list [GeV^4]", (Float64 ∘ EUval).(GeV, ρ_SM_list[dump_indices]))
        write(io, "m_PBH_list [gram]", Float64.(m_PBH_list[dump_indices] ./ NU.g))
        write(io, "n_PBH_list [cm^-3]", Float64.(n_PBH_list[dump_indices] .* (1e-2 * NU.m)^3))
        write(io, "Xₙ_list", Float64.(Xₙ_list[dump_indices]))
        write(io, "T_SM_list [GeV]", (Float64 ∘ EUval).(GeV, T_SM_list[dump_indices]))
        write(io, "Yₚ", Float64(Yₚ))
        if mesons_output_flag
            write(io, "n_pion_minus_list [cm^-3]", Float64.(n_pion_minus_list[dump_indices] .* (1e-2 * NU.m)^3))
            write(io, "n_pion_plus_list [cm^-3]", Float64.(n_pion_plus_list[dump_indices] .* (1e-2 * NU.m)^3))
            write(io, "n_kaon_minus_list [cm^-3]", Float64.(n_kaon_minus_list[dump_indices] .* (1e-2 * NU.m)^3))
            write(io, "n_kaon_zero_long_list [cm^-3]", Float64.(n_kaon_zero_long_list[dump_indices] .* (1e-2 * NU.m)^3))
        end
    end

    mesons_output_flag &&
        return a_list, t_list, ρ_SM_list, m_PBH_list, n_PBH_list, Xₙ_list, Yₚ,
               n_pion_minus_list, n_pion_plus_list, n_kaon_minus_list, n_kaon_zero_long_list

    return a_list, t_list, ρ_SM_list, m_PBH_list, n_PBH_list, Xₙ_list, Yₚ
end

@info "Loading `main-scan-parameters.jl`"
include("main-scan-parameters.jl")

Yₚ_mat = pmap(params -> (last ∘ magic)(params...), Iterators.product(m_PBH_f_list, β_PBH_list))
χ²_mat = ((Yₚ_mat .- .245) ./ .003).^2
jldopen(joinpath(output_data_directory, "scanning.jld2"), "w+"; compress=true) do io
    io["m_PBH_list [gram]"] = Float64.(m_PBH_f_list ./ NU.g)
    io["β_PBH_list"] = Float64.(β_PBH_list)
    io["Y_P_mat"] = Float64.(Yₚ_mat)
    io["chi2_mat"] = Float64.(χ²_mat)
end

@info "Scanned Yₚ and saved to $(joinpath(dump_data_directory, "scanning.jld2"))."
