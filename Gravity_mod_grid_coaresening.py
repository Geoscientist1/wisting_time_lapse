import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from skimage.measure import block_reduce
# from dsw_dpt_slices import GRID, POR, dz, fluid_type, typo
# from files_readme import celltopdpt, phi, cell_vol
from files_readme import sat27, sat29, sat31, sat32, sat33, sat35, sat39, sat42, sat45, sat49, sat53, sat57
# from files_readme import sgas_27, sgas_29, sgas_31, sgas_33, sgas_35, sgas_39, sgas_45, sgas_49, sgas_53, sgas_57
# from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
# from files_readme import A, B, phi, cel_cen_dpt, dx, dy, dz, dx_accum, dy_accum

rho27_27 = np.loadtxt('rho27_27.txt', delimiter=' ')
rho27_29 = np.loadtxt('rho27_29.txt', delimiter=' ')
rho27_31 = np.loadtxt('rho27_31.txt', delimiter=' ')
rho27_33 = np.loadtxt('rho27_33.txt', delimiter=' ')
rho27_35 = np.loadtxt('rho27_35.txt', delimiter=' ')
rho27_39 = np.loadtxt('rho27_39.txt', delimiter=' ')
rho27_45 = np.loadtxt('rho27_45.txt', delimiter=' ')
rho27_49 = np.loadtxt('rho27_49.txt', delimiter=' ')
rho27_53 = np.loadtxt('rho27_53.txt', delimiter=' ')
rho27_57 = np.loadtxt('rho27_57.txt', delimiter=' ')
########################################################################################################################


# dx = np.loadtxt("dx_uniform.txt", delimiter=" ")
# dx_accum = np.loadtxt("dx_accum_uniform.txt", delimiter=" ")
# dy = np.loadtxt("dy_uniform.txt", delimiter=" ")
# dy_accum = np.loadtxt("dy_accum_uniform.txt", delimiter=" ")

rho_oil = 835  # Oil density
rho_gas = 0.65  # Gas density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
k_range = 101
k_range_new = 5  # (2, 5, 10)
mod_fk = 20  # (50, 20, 10)
rho0 = 2650
# dv = dx * dy * dz
inv_B = 80
inv_A = 79
B = 159
A = 157

# This code makes new coarser grids out of the fine and complex petrel simulation grids, to speed up
# the gravity forward modeling, to be used in gravity inversion
################################################################################
################################################################################

def uniform_horizontal_gridding():
    m = 0
    dx_m = np.zeros(k_range * A * B)
    dy_m = np.zeros(k_range * A * B)
    dz_m = np.zeros(k_range * A * B)
    dx_accum = np.zeros(k_range * A * B)
    dy_accum = np.zeros(k_range * A * B)
    for k in range(k_range):
        for j in range(B):
            for i in range(A):
                if dx[m] == 0 or dy[m] == 0 or dz[m] == 0:
                    m += 1
                    continue
                else:
                    dx_m[m] = 70
                    dy_m[m] = 70
                    dz_m[m] = 2
                    m += 1
    i = 0
    m = 0
    temp = 0
    count_dx = 0
    count_dy = 0
    for k in range(k_range):
        for j in range(B):
            if dy[m] == 0 and i == 0:
                count_dy += 70
            else:
                dy_accum[m] = dy_m[m] + count_dy
                temp = dy_accum[m]
                count_dy += 70
            for i in range(A):
                if dx_m[m] == 0:
                    count_dx += 70
                    m += 1
                else:
                    dx_accum[m] = dx_m[m] + count_dx
                    dy_accum[m] = temp
                    count_dx += 70
                    m += 1
        return dx_m, dy_m, dz_m, dx_accum, dy_accum


# def density_modelling(sat_1, sat_2, soil_1, soil_2, sgas_1, sgas_2):
#     # rho_1 = ((1 - soil_1) * rho_water) + ((1 - sat_1) * rho_oil)
#     # rho_1 = (sat_1 * rho_water) + (soil_1 * rho_oil) + (sgas_1 * rho_gas)
#     # rho_2 = (sat_2 * rho_water) + (soil_2 * rho_oil) + (sgas_2 * rho_gas)
#     # rho_2 = ((1 - soil_2) * rho_water) + ((1 - sat_2) * rho_oil)
#
#     rho_1 = (sat_1 * rho_water) + ((1 - sat_1) * rho_oil)
#     rho_2 = (sat_2 * rho_water) + ((1 - sat_2) * rho_oil)
#
#     delta_rho = rho_2 - rho_1
#     return delta_rho
#
#
# rho27_27 = density_modelling(sat27, sat27, soil_27, soil_27, sgas_27, sgas_27)
# np.savetxt('rho27_27.txt', rho27_27, delimiter=' ', header='')
# rho27_29 = density_modelling(sat27, sat29, soil_27, soil_29, sgas_27, sgas_29)
# np.savetxt('rho27_29.txt', rho27_29, delimiter=' ', header='')
# rho27_31 = density_modelling(sat27, sat31, soil_27, soil_31, sgas_27, sgas_31)
# np.savetxt('rho27_31.txt', rho27_31, delimiter=' ', header='')
# rho27_33 = density_modelling(sat27, sat33, soil_27, soil_33, sgas_27, sgas_33)
# np.savetxt('rho27_33.txt', rho27_33, delimiter=' ', header='')
# rho27_35 = density_modelling(sat27, sat35, soil_27, soil_35, sgas_27, sgas_35)
# np.savetxt('rho27_35.txt', rho27_35, delimiter=' ', header='')
# rho27_39 = density_modelling(sat27, sat39, soil_27, soil_39, sgas_27, sgas_39)
# np.savetxt('rho27_39.txt', rho27_39, delimiter=' ', header='')
# rho27_45 = density_modelling(sat27, sat45, soil_27, soil_45, sgas_27, sgas_45)
# np.savetxt('rho27_45.txt', rho27_45, delimiter=' ', header='')
# rho27_49 = density_modelling(sat27, sat49, soil_27, soil_49, sgas_27, sgas_49)
# np.savetxt('rho27_49.txt', rho27_49, delimiter=' ', header='')
# rho27_53 = density_modelling(sat27, sat53, soil_27, soil_53, sgas_27, sgas_53)
# np.savetxt('rho27_53.txt', rho27_53, delimiter=' ', header='')
# rho27_57 = density_modelling(sat27, sat57, soil_27, soil_57, sgas_27, sgas_57)
# np.savetxt('rho27_57.txt', rho27_57, delimiter=' ', header='')


def mod_grid_coarsening_dz(dt):
    dt_coarse = np.zeros((k_range_new, B, A))
    dt = np.reshape(dt, (k_range, B, A))
    for i in range(A):
        for j in range(B):
            s = 0
            temp = 0
            for k in range(k_range):
                # if dt[k, j, i] != 0:
                #     X = dt[k, j, i]
                temp += dt[k, j, i]
                if k != 0 and mod(k, mod_fk) == 0:
                    dt_coarse[s, j, i] = temp
                    s += 1
                    temp = 0
    result = dt_coarse.flatten()
    return result


########################################################################################################################
# This is only to use when the row 2 in slide 311 (presentasjon 1, figures 5, 6, 7 and 8) method is to be applied
# def mod_grid_coarsening_dz(dt):
#     dt_coarse = np.zeros((k_range_new, B, A))
#     dt = np.reshape(dt, (k_range, B, A))
#     for i in range(A):
#         for j in range(B):
#             s = 0
#             temp = 0
#             for k in range(k_range):
#                 if dt[k, j, i] == 0:
#                     summen = 0
#                     count = 0
#                     for n in range(k_range):
#                         summen += dt[n, j, i]
#                         if dt[n, j, i] != 0:
#                             count += 1
#                     if count != 0:
#                         dt[k, j, i] = summen / count
#                 temp += dt[k, j, i]
#                 if k != 0 and mod(k, mod_fk) == 0:
#                     dt_coarse[s, j, i] = temp
#                     s += 1
#                     temp = 0
#     result = dt_coarse.flatten()
#     return result
########################################################################################################################

def mod_grid_coarsening(dt):
    dt_coarse = np.zeros((k_range_new, B, A))
    dt = np.reshape(dt, (k_range, B, A))
    for i in range(A):
        for j in range(B):
            s = 0
            fk = 0
            temp = 0
            for k in range(k_range):
                if dt[k, j, i] != 0:
                    fk += 1
                    X = dt[k, j, i]
                    temp += dt[k, j, i]
                if k != 0 and mod(k, mod_fk) == 0 and fk != 0:
                    Y = dt[k, j, i]
                    dt_coarse[s, j, i] = temp / fk
                    s += 1
                    temp = 0
                    fk = 0
    result = dt_coarse.flatten()
    for i in range(len(result)):
        if result[i] == 0.0:
            result[i] = 1.0
    return result


########################################################################################################################
# This is only to use when the row 2 in slide 311 (presentasjon 1, figures 5, 6, 7 and 8) method is to be applied
# def mod_grid_coarsening(dt):
#     dt_coarse = np.zeros((k_range_new, B, A))
#     dt = np.reshape(dt, (k_range, B, A))
#     for i in range(A):
#         for j in range(B):
#             s = 0
#             fk = 0
#             temp = 0
#             for k in range(k_range):
#                 if dt[k, j, i] == 0:
#                     summen = 0
#                     count = 0
#                     for n in range(k_range):
#                         summen += dt[n, j, i]
#                         if dt[n, j, i] != 0:
#                             count += 1
#                     if count != 0:
#                         dt[k, j, i] = summen / count
#                         X = dt[k, j, i]
#                 if dt[k, j, i] != 0:
#                     fk += 1
#                 temp += dt[k, j, i]
#                 if k != 0 and mod(k, mod_fk) == 0 and fk != 0:
#                     Y = dt[k, j, i]
#                     dt_coarse[s, j, i] = temp / fk
#                     s += 1
#                     temp = 0
#                     fk = 0
#     result = dt_coarse.flatten()  # Remember to change back to "result = dt_coarse.flatten()" when you are doing the
#     # averaging again
#     return result
########################################################################################################################


def mod_grid_coarsening_rho(dt):  # This is only to use when the adding and substracting rho0 method is to be applied
    dt += rho0
    dt_coarse = np.zeros((k_range_new, B, A))
    dt = np.reshape(dt, (k_range, B, A))
    for i in range(A):
        for j in range(B):
            s = 0
            fk = 0  # not needed
            temp = 0
            for k in range(k_range):
                if dt[k, j, i] != 0:
                    fk += 1  # This will never run since all the values in the different density arrays are non-zero.
                temp += dt[k, j, i]
                if k != 0 and mod(k, mod_fk) == 0 and fk != 0:  # you can also remove fk here
                    Y = dt[k, j, i]
                    dt_coarse[s, j, i] = temp / mod_fk  # not fk anymore here!
                    s += 1
                    temp = 0
                    fk = 0  # not needed
    result = dt_coarse.flatten()
    return result


def rho0_array_generator():  # This is only to use when the adding and substracting rho0 method is to be applied
    dt_coarse = np.zeros((k_range_new, B, A))
    for i in range(A):
        for j in range(B):
            for k in range(k_range_new):
                dt_coarse[k, j, i] = rho0
    result = dt_coarse.flatten()
    return result


########################################################################################################################


########################################################################################################################


def inv_grid_generator(dt):  # A FUNCTION FOR GENERATING A 79*80*5 INVERSION GRID FOR FASTER INVERSION RUNS ON THE
    # CLUSTER

    dt_coarse = np.zeros((k_range_new, inv_B, A))
    dt = np.reshape(dt, (k_range_new, B, A))
    modx = 2
    for i in range(A):
        s = 0
        fk = 0
        temp = 0
        for j in range(B):
            for k in range(k_range_new):
                if dt[k, j, i] != 1.0:
                    fk += 1
                    temp += dt[k, j, i]
                if j != 0 and mod(j, modx) == 0 and fk != 0:
                    dt_coarse[k, s, i] = temp / fk
                    s += 1
                    temp = 0
                    fk = 0

    dt_coarse_new = np.zeros((k_range_new, inv_B, inv_A))
    s = 0
    fk = 0
    temp = 0
    for i in range(A):
        for j in range(inv_B):
            for k in range(k_range_new):
                if dt[k, j, i] != 1.0:
                    fk += 1
                    temp += dt_coarse[k, j, i]
                if i != 0 and mod(i, modx) == 0 and fk != 0:
                    dt_coarse_new[k, j, s] = temp / fk
                    s += 1
                    temp = 0
                    fk = 0
                if i == 156:
                    s = 0
                    fk = 0
                    temp = 0
    result = dt_coarse_new
    return result


########################################################################################################################

# dx_uni, dy_uni, dz_uni, dx_accum_uni, dy_accum_uni = uniform_horizontal_gridding()

# np.savetxt('dx_uniform.txt', dx_uni, delimiter=' ', header='Uniform DX values')
# np.savetxt('dy_uniform.txt', dy_uni, delimiter=' ', header='Uniform DY values')
# np.savetxt('dz_uniform.txt', dz_uni, delimiter=' ', header='Uniform DZ values')
# np.savetxt('dx_accum_uniform.txt', dx_accum_uni, delimiter=' ', header='Uniform DX_accum values')
# np.savetxt('dy_accum_uniform.txt', dy_accum_uni, delimiter=' ', header='Uniform Dy_accum values')
########################################################################################################################
# This is to choose if the grid is uniform (dx = dy = 70m, and dz = 2m)
# dx_coarse = mod_grid_coarsening(dx_uni)
# dy_coarse = mod_grid_coarsening(dy_uni)
# dz_coarse = mod_grid_coarsening(dz_uni)
# dx_accum_coarse = mod_grid_coarsening(dx_accum)
# dy_accum_coarse = mod_grid_coarsening(dy_accum)
########################################################################################################################
# This is to choose if the grid is non-uniform (Petrel grid values for dx, dy and dz)
# dx_coarse = mod_grid_coarsening(dx)
# dy_coarse = mod_grid_coarsening(dy)
# dz_coarse = mod_grid_coarsening_dz(dz)
# dx_accum_coarse = mod_grid_coarsening(dx_accum)
# dy_accum_coarse = mod_grid_coarsening(dy_accum)
#######################################################################################################################
# dv_coarse = mod_grid_coarsening(dv)
# cel_cen_dpt_coarse = mod_grid_coarsening(cel_cen_dpt)
# phi_coarse = mod_grid_coarsening(phi)
# rho0_coarse = rho0_array_generator()
#######################################################################################################################
# sat27_coarse = mod_grid_coarsening(sat27)
sat32_coarse = mod_grid_coarsening(sat32)
sat42_coarse = mod_grid_coarsening(sat42)
# sat57_coarse = mod_grid_coarsening(sat57)
# test = np.reshape(sat57_coarse, (5, 159, 157))
# soil_27_coarse = mod_grid_coarsening(soil_27)
# soil_57_coarse = mod_grid_coarsening(soil_57)
#######################################################################################################################
# This is only to use when the adding and substracting rho0 method is to be applied
# rho27_27_coarse = mod_grid_coarsening_rho(rho27_27)
# rho27_29_coarse = mod_grid_coarsening_rho(rho27_29)
# rho27_31_coarse = mod_grid_coarsening_rho(rho27_31)
# rho27_33_coarse = mod_grid_coarsening_rho(rho27_33)
# rho27_35_coarse = mod_grid_coarsening_rho(rho27_35)
# rho27_39_coarse = mod_grid_coarsening_rho(rho27_39)
# rho27_45_coarse = mod_grid_coarsening_rho(rho27_45)
# rho27_49_coarse = mod_grid_coarsening_rho(rho27_49)
# rho27_53_coarse = mod_grid_coarsening_rho(rho27_53)
# rho27_57_coarse = mod_grid_coarsening_rho(rho27_57)
#######################################################################################################################
# rho27_27_coarse = mod_grid_coarsening(rho27_27)
# rho27_29_coarse = mod_grid_coarsening(rho27_29)
# rho27_31_coarse = mod_grid_coarsening(rho27_31)
# rho27_33_coarse = mod_grid_coarsening(rho27_33)
# rho27_35_coarse = mod_grid_coarsening(rho27_35)
# rho27_39_coarse = mod_grid_coarsening(rho27_39)
# rho27_45_coarse = mod_grid_coarsening(rho27_45)
# rho27_49_coarse = mod_grid_coarsening(rho27_49)
# rho27_53_coarse = mod_grid_coarsening(rho27_53)
# rho27_57_coarse = mod_grid_coarsening(rho27_57)

# np.savetxt('dx_coarse.txt', dx_coarse, delimiter=' ', header='')
# np.savetxt('dy_coarse.txt', dy_coarse, delimiter=' ', header='')
# np.savetxt('dz_coarse.txt', dz_coarse, delimiter=' ', header='')
# np.savetxt('dv_coarse.txt', dv_coarse, delimiter=' ', header='')
# np.savetxt('dx_accum_coarse.txt', dx_accum_coarse, delimiter=' ', header='')
# np.savetxt('dy_accum_coarse.txt', dy_accum_coarse, delimiter=' ', header='')
# np.savetxt('cel_cen_dpt_coarse.txt', cel_cen_dpt_coarse, delimiter=' ', header='')
# np.savetxt('phi_coarse.txt', phi_coarse, delimiter=' ', header='')

# np.savetxt('sat27_coarse.txt', sat27_coarse, delimiter=' ', header='')
np.savetxt('sat32_coarse.txt', sat32_coarse, delimiter=' ', header='')
np.savetxt('sat42_coarse.txt', sat42_coarse, delimiter=' ', header='')
# np.savetxt('sat57_coarse.txt', sat57_coarse, delimiter=' ', header='')
# np.savetxt('soil_27_coarse.txt', soil_27_coarse, delimiter=' ', header='')
# np.savetxt('soil_57_coarse.txt', soil_57_coarse, delimiter=' ', header='')

# np.savetxt('rho0_coarse.txt', rho0_coarse, delimiter=' ', header='')
# np.savetxt('rho27_27_coarse.txt', rho27_27_coarse, delimiter=' ', header='')
# np.savetxt('rho27_29_coarse.txt', rho27_29_coarse, delimiter=' ', header='')
# np.savetxt('rho27_31_coarse.txt', rho27_31_coarse, delimiter=' ', header='')
# np.savetxt('rho27_33_coarse.txt', rho27_33_coarse, delimiter=' ', header='')
# np.savetxt('rho27_35_coarse.txt', rho27_35_coarse, delimiter=' ', header='')
# np.savetxt('rho27_39_coarse.txt', rho27_39_coarse, delimiter=' ', header='')
# np.savetxt('rho27_45_coarse.txt', rho27_45_coarse, delimiter=' ', header='')
# np.savetxt('rho27_49_coarse.txt', rho27_49_coarse, delimiter=' ', header='')
# np.savetxt('rho27_53_coarse.txt', rho27_53_coarse, delimiter=' ', header='')
# np.savetxt('rho27_57_coarse.txt', rho27_57_coarse, delimiter=' ', header='')
