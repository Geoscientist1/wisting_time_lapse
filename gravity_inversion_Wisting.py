import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from Grav_Modeling_coarse_grid import grav_mod_coarse, grav_mod_coarse_real, grav_mod_coarse_final, \
    grav_mod_coarse_initial, grav_mod_coarse_inversion
from Grav_Modeling_coarse_grid import XSCALE, YSCALE, sat27, sat57, phi, cel_cen_dpt, dx_accum, dy_accum, dx, dy, dz
import random as rn
from numpy.linalg import norm

# from dsw_dpt_slices import GRID, POR, dz, k_range, fluid_type, typo
# from files_readme import celltopdpt, phi, cell_vol
# from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
# from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
# from files_readme import A, B, phi, cel_cen_dpt, dx, dy, dz, dx_accum, dy_accum
########################################################################################################################
rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
R_V = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
g_constant = 6.67 * 1e-11  # Gravity constant
int_min = 0
int_max = 0.5
rsv = 0.5
bck = 0.1  # Background water saturation time-lapse change
x_range = 157  # This is the # of grid cells in the x-direction
y_range = 159  # This is the # of grid cells in the y-direction
z_range = 5  # This is the # of grid cells in the z-direction
x_range_data = 55  # This is the # of grid data points  in the x-direction
y_range_data = 57  # This is the # of grid data points  in the y-direction
x_range_inv_grid = 117  # This is the # of grid cells in the inversion grid in the x-direction
y_range_inv_grid = 119  # This is the # of grid cells in the inversion grid in the y-direction
x_max = 13213.48
y_max = 14173.40
z_max = -1200
N = x_range_data * y_range_data  # This is the # number of data points/measurement stations
M = x_range * y_range * z_range  # This is the # number of data model parameters
S = x_range_inv_grid * y_range_inv_grid * z_range  # This is the # number of data model parameters in the inversion grid
# delta_SW_starting_model = "random between 0 and 0.5"  # too conservative guess
delta_SW_starting_model = 0.5  # too conservative guess
mu_r = [9000]  # Regularization parameter
fig_num = 183


########################################################################################################################
# This is a script that computes the starting rho0 density or model vector to be used in the inversion scheme
def dws_start_model(dsw):
    dsw0 = (np.ones(x_range * y_range * z_range)) * dsw
    return dsw0


def dws_start_model_random():
    dsw0 = (np.ones(x_range * y_range * z_range))
    for i in range(len(dsw0)):
        dsw0[i] = rn.uniform(int_min, int_max)
    return dsw0


########################################################################################################################
# This is a script that computes the necessary quantities to run gravity inversion inshAllah
def real_density_distribution(sat_1, sat_2):
    rho_1 = (sat_1 * rho_water) + ((1 - sat_1) * rho_oil)
    rho_2 = (sat_2 * rho_water) + ((1 - sat_2) * rho_oil)
    delta_rho = rho_2 - rho_1
    rho = delta_rho * phi
    return rho


def gravity_data_misfit(observed_data, modeled_data):
    v = np.zeros((D, D))
    for i in range(D):
        for j in range(D):
            if i == j:
                v[j, i] = 1 / (0.001 + (0.03 * observed_data[j, i]))
    misfit = (1 / 2) * norm(v * (modeled_data - observed_data))
    return misfit


def gravity_kernel(cell_dpt_input, dx_accum_input, dy_accum_input, length, x_range_d, y_range_d, x_scale_input,
                   y_scale_input):
    G_kernel = np.zeros((x_range_d * y_range_d, length))
    a = (6.67 * 1e-11) * (cell_dpt_input - 400)
    n = 0
    print("Kernel computation...")
    for y_step in range(0, y_range_d):
        for x_step in range(0, x_range_d):
            m = 0
            # print("Kernel computation")
            # print("y", y_step)
            # print("x", x_step)
            for k in range(length):
                # r = np.sqrt((((x_step * x_scale_input) + (x_scale_input / 2)) - dx_accum_input[m]) ** 2 + (
                #         ((y_step * y_scale_input) + (y_scale_input / 2)) - dy_accum_input[m]) ** 2)
                r = np.sqrt((((x_step * x_scale_input) + x_scale_input) - dx_accum_input[m]) ** 2 + (
                        ((y_step * y_scale_input) + y_scale_input) - dy_accum_input[m]) ** 2)
                d = pow((pow(r, 2) + pow(cell_dpt_input[m] - 400, 2)), 1.5)
                G_kernel[n, m] = ((a[m]) / d)
                m += 1
            n += 1
    return G_kernel


def data_weighting_matrix(flr, percent, observed_data):
    v = np.zeros((N, N))
    observed_data = observed_data.flatten()
    for j in range(N):
        for i in range(N):
            if i == j:
                # v[j, i] = 1 / (0.0001 * 1e-8 + (0.02 * observed_data[i]) * 1e-8)
                v[j, i] = 1 / (flr * 1e-8 + (percent * observed_data[i]) * 1e-8)
    return v


def depth_weighting_matrix(dimension, z, z0, beta):
    w_model = np.zeros((dimension, dimension))
    for i in range(dimension):
        for j in range(dimension):
            if i == j:
                w_model[j, i] = pow((z[i] + z0), -beta / 2)
    return w_model


def mod_parm_mapping(dswat_input):
    dsw_diff_real = dswat_input
    p = []
    por = []
    cell_depth = []
    dx_accumulated = []
    dy_accumulated = []
    dxx = []
    dyy = []
    dzz = []
    dt_diff = np.subtract(sat57, sat27)
    for i in range(M):
        if sat27[i] and sat57[i] != 1:
        # if sat27[i] and sat57[i] and phi[i] and cel_cen_dpt[i] != 1:
            p.append(dsw_diff_real[i])
            por.append(phi[i])
            cell_depth.append(cel_cen_dpt[i])
            dx_accumulated.append(dx_accum[i])
            dy_accumulated.append(dy_accum[i])
            dxx.append(dx[i])
            dyy.append(dy[i])
            dzz.append(dz[i])

    p = np.array(p)
    por = np.array(por)
    cell_depth = np.array(cell_depth)
    dx_accumulated = np.array(dx_accumulated)
    dy_accumulated = np.array(dy_accumulated)
    dxx = np.array(dxx)
    dyy = np.array(dyy)
    dzz = np.array(dzz)
    T = len(p)

    mapping_mtr = np.zeros((M, T))
    j = 0
    for i in range(M):
        if sat27[i] and sat57[i] != 1:
            mapping_mtr[i, j] = 1
            j += 1
    return p, por, cell_depth, dx_accumulated, dy_accumulated, dxx, dyy, dzz, T, mapping_mtr


########################################################################################################################
#         Conjugate Gradient (CG) Inversion Algorithm

def conjugate_gradient_inversion(dsw0, mu, g, v, w_weight, d_observed, d_initial, value, porr,
                                 cell_depthh, dx_accumulatedd, dy_accumulatedd, dxxx, dyyy, dzzz, map_mat):
    information = "None"
    d_initial = d_initial.flatten() * 1e-8
    d_observed = d_observed.flatten() * 1e-8
    misfit_r = (1 / 2) * norm(v * (np.subtract(d_initial, d_observed)))
    model_norm_initial = (1 / 2) * norm(w_weight * dsw0)
    # model_norm_initial = (1 / 2) * norm(dsw0)
    initial_total_misfit = misfit_r + model_norm_initial
    old_misfit = 0
    g_transpose = g.transpose()  # Please look at slide 360
    v_transpose = v.transpose()  # Please look at slide 360
    w_transpose = w_weight.transpose()  # Please look at slide 360
    A = (((g_transpose @ v_transpose) @ v) @ g) + (mu * (w_transpose @ w_weight))
    # A = (((g_transpose @ v_transpose) @ v) @ g) + mu
    b = (((g_transpose @ v_transpose) @ v) @ d_observed)
    r0 = (A @ dsw0) - b
    p0 = -r0
    dsw_k = dsw0
    rk = r0
    pk = p0

    count = 0
    iteration = [0]
    misfit_iteration = [initial_total_misfit]
    print("********************************")
    print("This is the initial data misfit:", misfit_r)
    print("This is the initial model norm:", model_norm_initial)
    print("This is the initial total misfit:", initial_total_misfit)
    print("********************************")

    while rk.all() != 0:
        # if 2500 < misfit_r < 3500:  # Target misfit tolerance
        if 100 < misfit_r < 150:
            information = "Target misfit tolerance is achieved!"
            break

        alpha = (rk.transpose() @ rk) / ((pk.transpose() @ A) @ pk)
        rho_new = dsw_k + (alpha * pk)
        r_new = rk + (alpha * (A @ pk))
        beta_r = (r_new.transpose() @ r_new) / (rk.transpose() @ rk)
        p_new = -r_new + (beta_r * pk)

        dsw_k = rho_new
        rk = r_new
        pk = p_new

        dsw_k_mapped = map_mat @ dsw_k
        # por_mapped = map_mat @ porr
        # cell_depth_mapped = map_mat @ cell_depthh
        # dx_accum_mapped = map_mat @ dx_accumulatedd
        # dy_accum_mapped = map_mat @ dy_accumulatedd
        # dx_mapped = map_mat @ dxxx
        # dy_mapped = map_mat @ dyyy
        # dz_mapped = map_mat @ dzzz
        modeled_data = grav_mod_coarse_inversion(dsw_k_mapped, count)
        ################################################################################################################
        np.savetxt('recovered_model_temp_' + str(fig_num) + '.txt', dsw_k, delimiter=' ', header='')
        np.savetxt('recovered_data_temp_' + str(fig_num) + '.txt', modeled_data, delimiter=' ', header='')
        ################################################################################################################
        modeled_data = modeled_data.flatten() * 1e-8
        misfit_temp = (1 / 2) * norm(v * (np.subtract(modeled_data, d_observed)))
        # model_norm_temp = (1 / 2) * norm(w_weight * (dsw_k_mapped))
        model_norm_temp = (1 / 2) * norm(w_weight * (dsw_k))
        # model_norm_temp = (1 / 2) * norm((dsw_k))
        total_misfit = misfit_temp + model_norm_temp
        misfit_r = misfit_temp

        if old_misfit == misfit_r:
            information = "The old and new misfits are identical!"
            break
        old_misfit = misfit_r
        count += 1
        iteration.append(count)
        misfit_iteration.append(total_misfit)

        print("This is the value of the data misfit:", misfit_r)
        ################################################################################################################
        plt.plot(iteration, misfit_iteration)
        plt.xlabel("Number of Iterations")
        plt.ylabel("Total Misfit ($\phi_m + \phi_d$)")
        ax = plt.gca()
        plt.savefig("inversion_iteration_" + str(fig_num))
        ################################################################################################################
    plt.plot(iteration, misfit_iteration, label="\u03bc= " + str(mu))
    plt.xlabel("Number of Iterations")
    plt.ylabel("Total Misfit ($\phi_m + \phi_d$)")
    plt.yscale("log")
    ax = plt.gca()
    legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    plt.legend(loc="upper right")
    # Put a nicer background color on the legend.
    legend.get_frame()
    plt.savefig("inversion_iteration_" + str(fig_num))
    ####################################################################################################################
    model_norm = (1 / 2) * norm(w_weight * (dsw_k))
    # model_norm = (1 / 2) * norm(dsw_k)
    print("The Conjugate Gradient Inversion is successfully done Alhamdulilah")
    print(information)
    return dsw_k, dsw0, misfit_r, model_norm, pk, iteration, misfit_iteration


def gravity_inversion_execution_Wisting():
    ####################################################################################################################

    #         FUNCTION EXECUTIONS

    ####################################################################################################################
    ####################################################################################################################
    #         STARTING MODEL
    ####################################################################################################################
    print(
        "*************************************************************************************************************")
    print(
        "WELCOME TO THE GRAVITY INVERSION SCRIPT, FROM HERE YOU WILL BE PROVIDED WITH dsw_0, dsw_observed, data_0 AND "
        "data_observed ")
    print(
        "*************************************************************************************************************")
    ####################################################################################################################
    # ####################################################################################################################
    # #         INITIAL MODEL
    # ####################################################################################################################
    dswat_start_model = dws_start_model(delta_SW_starting_model)
    # dswat_start_model = dws_start_model_random()
    print("This is the run number:", fig_num)
    print("This is the initial guess in water saturation:", delta_SW_starting_model)
    print("This is the regularization parameter:", mu_r)
    print("The FWM of the initial start model is running...")
    data_0 = grav_mod_coarse_initial(dswat_start_model)
    # np.savetxt('grav_data_57_start_model_Wisting_random.txt', data_0, delimiter=' ',
    #            header=' ')
    np.savetxt('grav_data_57_start_model_Wisting_' + str(fig_num) + '.txt', data_0, delimiter=' ',
               header=' ')
    # data_0 = np.loadtxt('grav_data_57_start_model_Wisting_0.7.txt', delimiter=' ')
    print("The FWM of the initial start model is done!")
    # ####################################################################################################################
    # #         OBSERVED MODEL
    # ####################################################################################################################
    dswat_observed = np.subtract(sat57, sat27)
    data_observed = np.loadtxt('grav_data_coarse_grid_coarse_dtpoint_grid.txt', delimiter=' ')
    # ####################################################################################################################
    # # Adding Gaussian noise to the synthetic observed data
    # ####################################################################################################################
    mud = 0
    sigma = 0.02
    max_value = 50
    min_value = 0
    noise = np.random.normal(mud, sigma, data_observed.shape)
    noise_data_observed = data_observed + noise
    # data_number = []
    # for i in range(x_range_data * y_range_data):
    #     data_number.append(i + 1)
    # d1 = data_observed.flatten()
    # d2 = noise_data_observed.flatten()
    # plt.plot(data_number[min_value:max_value], d1[min_value:max_value], "+", color='black', label="Accurate data")
    # plt.plot(data_number[min_value:max_value], d2[min_value:max_value], "o", fillstyle='none', color='black', label="Noisy data")
    # # y_error = d2
    # # plt.errorbar(data_number[0:10], d1[0:10],yerr=y_error[0:10], fmt='o', markersize=4, capsize=10)
    # plt.xlabel("data number")
    # plt.ylabel("Gravity data $(\u03BC$Gal)")
    # ax = plt.gca()
    # legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    # plt.legend(loc="upper right")
    # # Put a nicer background color on the legend.
    # legend.get_frame()
    # plt.show()
    ####################################################################################################################
    # plt.errorbar(data_number[min_value:max_value], d1[min_value:max_value], yerr=d2[min_value:max_value]-d1[min_value:max_value], fmt='o', markersize=4, capsize=10)
    # plt.xlabel("data number")
    # plt.ylabel("Gravity data $(\u03BC$Gal)")
    # ax = plt.gca()
    # legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    # plt.legend(loc="upper right")
    # # Put a nicer background color on the legend.
    # legend.get_frame()
    # plt.show()
    # np.savetxt('grav_data_coarse_grid_coarse_dtpoint_grid_noisy.txt', noise_data_observed, delimiter=' ',
    #            header='Noisy')
    # print("noise contamination is done!")
    ####################################################################################################################
    # Mapping model cells to inversion parameters
    ####################################################################################################################
    ####################################################################################################################
    new_dsw_diff, por, cell_depth, dx_accumulated, dy_accumulated, dxx, dyy, dzz, map_len, map_mat = \
        mod_parm_mapping(dswat_start_model)

    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################
    # Conjugate Gradient Inversion Analysis
    print("Conjugate Gradient Inversion Analysis starts now!")
    ####################################################################################################################
    flor = 0.0001
    percent = 0.02
    beta = 2
    z0 = 400
    # z0 = -375
    z = [600, 700, 800]
    # mu_r = 4.6
    NUMBER_OF_ITERATIONS = 1
    # mu_r = np.linspace(0.1, 1000, NUMBER_OF_ITERATIONS)
    # mu_r = [0.1, 1, 10, 100, 1000, 10000, 100000]
    misfit_d = np.zeros(NUMBER_OF_ITERATIONS)
    model_norm_d = np.zeros(NUMBER_OF_ITERATIONS)
    G_r = gravity_kernel(cell_depth, dx_accumulated, dy_accumulated, map_len, x_range_data, y_range_data, XSCALE,
                         YSCALE)
    V_r = data_weighting_matrix(flor, percent, noise_data_observed)
    W_r = depth_weighting_matrix(map_len, cell_depth, z0, beta)

    # Iterating through different and multiple values of the regularization parameter!

    recovered_model_sol, mu_sol, misfit_sol, model_norm_sol, iteration_number_sol, total_misfit_sol, solution_model = \
        None, None, None, None, None, None, None

    for i in range(NUMBER_OF_ITERATIONS):
        print("#######################################################################################################")
        print("Inversion iteration from scratch: ", i, " is running, for mu: ", mu_r[i])
        recovered_model, initial_model, misfit, model_norm_r, direction, iteration_number, misfit_iteration_result = \
            conjugate_gradient_inversion(new_dsw_diff, mu_r[i], G_r, V_r, W_r, noise_data_observed, data_0, map_len,
                                         por, cell_depth, dx_accumulated, dy_accumulated, dxx, dyy, dzz, map_mat)

        print("This is the initial model", dswat_start_model)
        print("This is the observed model", dswat_observed)
        print("This is the recovered model", recovered_model)
        print("This is the data misft:", misfit)
        print("This is the model norm:", model_norm_r)
        print("This is the regularization parameter:", mu_r[i])
        print("#######################################################################################################")
        misfit_d[i] = misfit
        model_norm_d[i] = model_norm_r
        total_misfit_d = misfit_iteration_result
        iteration_number_d = iteration_number

        # if 1000 < misfit_d[i] < 2000:
        # if 2500 < misfit_d[i] < 3500:
        # if 10000 < misfit_d[i] < 120000:
        recovered_model_sol = recovered_model
        mu_sol = mu_r[i]
        misfit_sol = misfit
        model_norm_sol = model_norm_r
        iteration_number_sol = iteration_number_d
        total_misfit_sol = total_misfit_d
        solution_model = map_mat @ recovered_model_sol
    #
    # ####################################################################################################################
    # # Plotting the data misfit as a function of the regularization parameter
    # plt.plot(mu_r, misfit_d)
    # plt.plot(mu_sol, misfit_sol, "*", color="tab:blue", markersize=7)
    # plt.axhline(y=9, color='orange', linestyle='-')
    # plt.xlabel("Regularization parameter " u"(\u03bc)")
    # plt.ylabel("Data misfit ($\phi_d$)")
    # # plt.ylim([0, 50])
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()
    # ####################################################################################################################
    # # Plotting the model norm as a function of the regularization parameter
    # plt.plot(mu_r, model_norm_d)
    # plt.plot(mu_sol, model_norm_sol, "*", color="tab:blue", markersize=7)
    # plt.xlabel("Regularization parameter " u"(\u03bc)")
    # plt.ylabel("Model norm ($\phi_m$)")
    # # plt.ylim([0, 50])
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()
    # ####################################################################################################################
    # Plotting the total misfit as a function of the iteration number
    # plt.plot(iteration_number_sol, total_misfit_sol)
    # plt.plot(mu_sol, model_norm_sol, "*", color="tab:blue", markersize=7)
    # plt.xlabel("Number of Iterations")
    # plt.ylabel("Total Misfit ($\phi_m + \phi_d$)")
    # plt.show()
    # ####################################################################################################################
    # ####################################################################################################################
    print(
        "*************************************************************************************************************")

    print("These are the solutions:", mu_sol, misfit_sol, model_norm_sol)
    print(
        "*************************************************************************************************************")
    # ####################################################################################################################
    solution_file_name = 'solution_details_' + str(fig_num) + '.txt'
    np.savetxt(solution_file_name, (mu_sol, misfit_sol, model_norm_sol), delimiter=' ',
               header='These are the solutions: mu_sol, data_misfit_sol, model_norm_sol ')
    print("Forward Modeling computation of the recovered model...")
    recovered_data = grav_mod_coarse_final(solution_model)

    return dswat_start_model, dswat_observed, data_observed, solution_model, recovered_data

########################################################################################################################
########################################################################################################################
#         Data Misfit
########################################################################################################################
# d_initial = np.loadtxt('grav_data_57_start_model_Wisting.txt', delimiter=' ')
# d_observed = np.loadtxt('grav_data_coarse_grid_coarse_dtpoint_grid_noisy.txt', delimiter=' ')
# flooor = 0.1
# percent = 0.02
# v = data_weighting_matrix(57, 55, flooor, percent, d_observed)
# data_misfit = (1 / 2) * norm(v * (np.subtract(d_initial, d_observed)))
# print("Alhamdulilah The misfit is:", data_misfit)
########################################################################################################################
########################################################################################################################
# Mapping simulation grid values to inversion parameters
########################################################################################################################
# dt_diff = np.subtract(sat57, sat27)
# p = []
# for i in range(M):
#     if dt_diff[i] != 0:
#         p.append(dt_diff[i])
# p = np.array(p)
#
# map_mtr = np.zeros((M, len(p)))
# j = 0
# for i in range(M):
#     if dt_diff[i] != 0:
#         map_mtr[i, j] = 1
#         j += 1
#
# print("Alhamdulilah")
#
# new_m = map_mtr @ p
#
# result = np.isin(dt_diff, new_m)
# print(result)

########################################################################################################################
########################################################################################################################
# Depth weighting analysis
########################################################################################################################
# alpha = np.linspace(0, 13213.48, 157)
# beta = np.linspace(0, 14173.40, 159)
# z = np.linspace(596.02, 1082.86, 5)
# x = np.linspace(0, 13213.48, 55)
# y = np.linspace(0, 14173.40, 57)
# h = 400
# w = np.zeros(5)
# kernel = np.zeros(M)
# z0 = -400
#
# dswat_start_model = dws_start_model(delta_SW_starting_model)
# dsw_diff = np.subtract(dswat_start_model, sat27)
# dsw_diff_real = dswat_start_model
# p = []
# por = []
# cell_depth = []
# dx_accumulated = []
# dy_accumulated = []
# dxx = []
# dyy = []
# dzz = []
# variable = -0.9
# for i in range(M):
#     if dsw_diff[i] != variable:
#         p.append(dsw_diff_real[i])
#         por.append(phi[i])
#         cell_depth.append(cel_cen_dpt[i])
#         dx_accumulated.append(dx_accum[i])
#         dy_accumulated.append(dy_accum[i])
#         dxx.append(dx[i])
#         dyy.append(dy[i])
#         dzz.append(dz[i])
#
# p = np.array(p)
# por = np.array(por)
# cell_depth = np.array(cell_depth)
# dx_accumulated = np.array(dx_accumulated)
# dy_accumulated = np.array(dy_accumulated)
# dxx = np.array(dxx)
# dyy = np.array(dyy)
# dzz = np.array(dzz)
# T = len(p)
# new_dsw_diff = p
#
# mapping_mtr = np.zeros((M, T))
# j = 0
# for i in range(M):
#     if dsw_diff[i] != variable:
#         mapping_mtr[i, j] = 1
#         j += 1
#
#
# for i in range(5):
#     w[i] = (1 / (z[i] + z0))
#
# m = 0
# for i in range(157):
#     for j in range(159):
#         for k in range(5):
#             r = np.sqrt(((x[0]) - alpha[i]) ** 2 + ((y[0]) - beta[j]) ** 2)
#             d = pow((pow(r, 2) + pow(z[k] - 400, 2)), 1.5)
#             kernel[m] = (z[k] - h) / d
#             m += 1
#             # kernel[i] = 1 / ((z[i] - h) ** 2)
#             m = 0
#
# G_kernel = np.zeros((57 * 55, T))
# n = 0
# for k in range(T):
#     r = np.sqrt((((x_step * x_scale_input) + (x_scale_input / 2)) - dx_accum_input[m]) ** 2 + (
#                         ((y_step * y_scale_input) + (y_scale_input / 2)) - dy_accum_input[m]) ** 2)
#     d = pow((pow(r, 2) + pow(cell_dpt_input[m] - 400, 2)), 1.5)
#     G_kernel[n, m] = ((a[m]) / d)
#     m += 1
# n += 1
#
# kernel = np.reshape(kernel, (157, 159, 5))
# print("Done")
# w = w ** 2
# plt.plot(w, z, label="Depth weighting")
# plt.plot(kernel[0][0], z, label="Kernel decay with depth")
# plt.plot(G_kernel[0][0], z, label="Kernel decay with depth")
# ax = plt.gca()
# ax.xaxis.tick_top()
# plt.gca().invert_yaxis()
# legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
# legend.get_frame()
# plt.xscale(value="log")
# plt.show()
# #
# plt.plot(w, z, label="Depth weighting")
# plt.plot(kernel[0], z, label="Kernel decay with depth")
# ax = plt.gca()
# plt.xscale(value="log")
# ax.xaxis.tick_top()
# plt.gca().invert_yaxis()
# legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
# plt.xlabel("Kernel value")
# plt.ylabel("Depth (m)")
# legend.get_frame()
# plt.show()
