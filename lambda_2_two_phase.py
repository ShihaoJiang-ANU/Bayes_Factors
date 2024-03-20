import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import scipy.stats as st


# F is the degree of melting in garnet field, F_2 in spinel field

def find_melt_composition(F, F2, Kp_grt, Kp_spl, K0_grt, K0_spl, c_mantle):
    c_residue = c_mantle * 1 / (1 - F) * (1 - Kp_grt * F / K0_grt) ** (1 / Kp_grt)
    F_prime = F2 / (1 - F)
    c_residue_two_phase = c_residue * 1 / (1 - F_prime) * (1 - Kp_spl * F_prime / K0_spl) ** (1 / Kp_spl)
    c_l = (c_mantle - c_residue_two_phase * (1 - F - F2)) / (F + F2)
    return c_l


PM_MS_1995 = {
    "La": 0.648, "Ce": 1.675, "Pr": 0.254, "Nd": 1.250, "Sm": 0.406, "Eu": 0.154,
    "Gd": 0.544, "Tb": 0.099, "Dy": 0.674, "Ho": 0.149, "Er": 0.438, "Tm": 0.068,
    "Yb": 0.441, "Lu": 0.0675}

c_mantle = np.array([0.648, 1.675, 1.250, 0.406, 0.154, 0.544, 0.674, 0.438, 0.441, 4.30])
# c_mantle = np.array([0.192, 0.55, 0.581, 0.239, 0.096, 0.358, 0.505, 0.348, 0.365, 3.328])

CI_REE = np.array([0.2347, 0.6032, 0.4524, 0.1471, 0.056, 0.1966, 0.2427, 0.1589, 0.1625, 1.56])

# Ol, Opx, Cpx, Spl, Gt from Michael J. Walter, 1998
mode_source_2GPa = [0.62, 0.281, 0.087, 0.011, 0.001]
mode_source_3GPa = [0.531, 0.177, 0.273, 0, 0.019]
mode_source_4GPa = [0.536, 0.056, 0.279, 0, 0.109]
mode_source_5GPa = [0.548, 0.000, 0.300, 0, 0.152]

ol_liq = [0.0000023, 0.0000073, 0.0000583, 0.0002905, 0.0005461,
          0.0009880, 0.0028912, 0.0066264, 0.0121092, 0.0039069]
opx_liq = [0.00074, 0.00159, 0.00600, 0.01579, 0.02266, 0.03145,
           0.05489, 0.08084, 0.10358, 0.06350]
cpx_liq = [0.055, 0.0876, 0.1878, 0.3083, 0.3638, 0.4169, 0.5034,
           0.5437, 0.5453, 0.5219]
sp_liq = [0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0006, 0.0015,
          0.003, 0.0045, 0.0045]
gt_liq = [0.001, 0.008, 0.057, 0.217, 0.45, 0.9, 2, 3.5, 7, 7]

k_data = np.array([ol_liq, opx_liq, cpx_liq, sp_liq, gt_liq])

#  sp + opx = gt + olv
#  MgAl2O4 (Spinel) + 2Mg2Si2O6 --> Mg3Al2Si3O12 (Garnet) +Mg2SiO4
#  142.265 + 200.78*2 = 403.1325 + 140.6925

reaction_with_cpx = np.array([-0.22, 0.38, 0.71, 0.13, 0])
reaction_with_gt = np.array([0.07, -0.16, 0.68, 0, 0.25])  # 3 GPa
# reaction_with_gt = np.array([0.10, -0.12, 0.58, 0, 0.25])  # 4 GPa
# reaction_with_gt = np.array([0.12, -0.10, 0.55, 0, 0.24])  # 6 GPa
# reaction_with_gt = np.array([0.26, 0, 0.50, 0, 0.24])  # 7 GPa

orth_value = np.array([
    [1, 0.1052, 4.913, 0.243, 8.29],
    [1, 0.0882, 2.033, -0.031, -7.16],
    [1, 0.0542, -1.994, -0.208, -1.99],
    [1, 0.0242, -3.627, -0.108, 8.64],
    [1, 0.0112, -3.776, -0.032, 9.67],
    [1, -0.0018, -3.587, 0.043, 7.91],
    [1, -0.0278, -2.194, 0.142, -1.45],
    [1, -0.0508, 0.165, 0.101, -7.63],
    [1, -0.0698, 2.912, -0.078, -0.16]
])
orth_value[:, 2] = orth_value[:, 2]*1e-3
orth_value[:, 3] = orth_value[:, 3]*1e-3
orth_value[:, 4] = orth_value[:, 4]*1e-6


def calculate_lambdas(F_gt, F_sp):
    mode_source = np.array(mode_source_4GPa)

    mode_gt_end = (mode_source - reaction_with_gt * F_gt) / np.sum(mode_source - reaction_with_gt * F_gt)
    # The following is incorrect
    # mode_gt_end = (mode_source - reaction_with_gt * F_gt) / (1 - F_gt)

    transition_reaction = np.array([-140.6925, 200.78 * 2, 0, 142.265, -403.1325])
    mode_sp_start = mode_gt_end + transition_reaction / 403.1325 * mode_gt_end[-1]
    # the following is incorrect
    # mode_sp_start = mode_gt_end + transition_reaction / 403.1325 * mode_gt_end

    k0_gt = np.dot(k_data.T, mode_source.T)
    kp_gt = np.dot(k_data.T, reaction_with_gt.T)
    k0_sp = np.dot(k_data.T, mode_sp_start.T)
    kp_sp = np.dot(k_data.T, reaction_with_cpx.T)

    melt_composition = find_melt_composition(F_gt, F_sp, kp_gt, kp_sp, k0_gt, k0_sp, c_mantle)
    log_value = np.log10(melt_composition/CI_REE)[0: 9]

    A = orth_value
    B = log_value.T

    # find x so that Ax = B
    x, residuals, rank, singular_values = np.linalg.lstsq(A, B, rcond=None)
    return x[2]


matrix = [[0 for _ in range(10)] for _ in range(10)]

for i in range(10):
    for j in range(10):
        matrix[i][j] = calculate_lambdas(i * 0.01, j * 0.01)

print(np.array(matrix))

