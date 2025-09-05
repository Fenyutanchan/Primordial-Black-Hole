# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

sqrt2 = sqrt(2)
ζ₃ = zeta(3)

EU = GeV
NU = NaturalUnit(EU)
BITS = 100

G_F = GeV(1.1663788e-5, -2) # Fermi constant
M_W = GeV(80.4335)
M_Z = GeV(91.1880)
m_e = MeV(0.51099895000)
m_μ = MeV(105.6583755)
Λ_QCD = MeV(200)
m_p = MeV(938.27208816)
m_n = MeV(939.5654205)
τ_n = 878.4 * NU.s
m_c = GeV(1.2730)
m_τ = MeV(1776.93)
m_b = GeV(4.183)

m_π = MeV(139.57039)
τ_π = 2.6033e-8 * NU.s
m_π0 = MeV(134.9768)
τ_π0 = 8.43e-17 * NU.s

m_K = MeV(493.677)
m_K0 = MeV(497.611)

f_π = MeV(130.2)
f_K = MeV(155.7)

V_ud = 0.97367
V_us = 0.22431

λ(x, y, z) = x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x

Q = m_n - m_p
α_EM = inv(137.035999084)
η_B = 6.1e-10
cubic_meter_per_second = NU.m^3 / NU.s
T_nuc = MeV(0.078185270688319950)
# T_nuc = MeV(0.07612105457906148)

Γ_pion_minus = Γ_pion_plus = inv(τ_π) # 3.8e7 / NU.s
Γ_kaon_minus = 8.1e7 / NU.s
Γ_kaon_zero_long = 2.0e7 / NU.s

function Fc(T::EnergyUnit, mₕ::EnergyUnit)
    vₑ = sqrt(T / mₕ) + sqrt(T / m_p)
    x = 2 * π * α_EM / vₑ
    return x / (1 - exp(-x))
end
Fc_pion_minus(T) = Fc(T, m_π)
Fc_kaon_minus(T) = Fc(T, m_K)

σv_pion_minus_pton(T) = 4.3e-23 * Fc_pion_minus(T) * cubic_meter_per_second
σv_pion_plus_ntop(T) = (4.3e-23 / 0.9) * cubic_meter_per_second
σv_kaon_minus_pton(T) = 9.6e-22 * Fc_kaon_minus(T) * cubic_meter_per_second
σv_kaon_minus_ntop(T) = (9.6e-22 / 2.46) * cubic_meter_per_second
σv_kaon_zero_long_pton(T) = 1.95e-22 * cubic_meter_per_second
σv_kaon_zero_long_ntop(T) = (1.95e-22 / 0.41) * cubic_meter_per_second
