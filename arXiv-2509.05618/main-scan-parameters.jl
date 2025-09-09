# Copyright (c) 2025 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

m_PBH_min_in_gram = 6e8
m_PBH_max_in_gram = 1e10
number_of_m_PBH = 100

β_PBH_min = 1e-20
β_PBH_max = 1e-17
number_of_β_PBH = 100

@info """
Scanning parameters:
m_PBH (min): $(m_PBH_min_in_gram) gram
m_PBH (max): $(m_PBH_max_in_gram) gram
number of m_PBH: $(number_of_m_PBH)

β_PBH (min): $(β_PBH_min)
β_PBH (max): $(β_PBH_max)
number of β_PBH: $(number_of_β_PBH)
"""

m_PBH_f_list = geomspace(ArbFloat(m_PBH_min_in_gram; bits=BITS), ArbFloat(m_PBH_max_in_gram; bits=BITS), number_of_m_PBH) .* NU.g
β_PBH_list = geomspace(ArbFloat(β_PBH_min; bits=BITS), ArbFloat(β_PBH_max; bits=BITS), number_of_β_PBH)
