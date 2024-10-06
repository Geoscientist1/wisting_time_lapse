import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from Grav_Modeling_coarse_grid import grav_mod_coarse, grav_mod_coarse_real, grav_mod_coarse_final, \
    grav_mod_coarse_initial, grav_mod_coarse_inversion
# from Grav_Modeling_coarse_grid import x_range_data, y_range_data, XSCALE, YSCALE, gv_kernel, data_observed, timestep,\
#     sat27, sat57, phi, cel_cen_dpt, dpt_cell, dx, dy, dz
from Grav_Modeling_coarse_grid import x_range_data, y_range_data, XSCALE, YSCALE, gv_kernel, timestep, data_observed, \
    sat27, sat32, sat57, cel_cen_dpt, dt_grid_resolution
from Grav_Gauss_FFT_FWM import gauss_fft_fwm
import random as rn
from numpy.linalg import norm
import time

# from dsw_dpt_slices import GRID, POR, dz, k_range, fluid_type, typo
# from files_readme import celltopdpt, phi, cell_vol
# from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
# from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
# from files_readme import A, B, phi, dpt, dx, dy, dz, dx_accum, dy_accum
########################################################################################################################
dpt = cel_cen_dpt
fig_num = 171
delta_SW_starting_model = 0.5  # too optimistic guess
int_min = 0.1
int_max = 0.7
z0 = -1300
mu_r = [1, 5, 10, 21, 50, 100, 200, 400, 600, 800, 1000, 2000, 4000, 7000, 10000, 50000, 100000]  # Regularization parameter
NUMBER_OF_ITERATIONS = 17
# -1500, -1300, -1100, -1000, -900, -800, -700, -600,-500, -300, -100, 0,  100, 200, 400, 600, 700
# 170,    171,   172,    173,  174,  175,   176,  177, 178   179   180  181 182  183  184  185  186

# -1300, -700,-900, -1100. -1000, -820,-690,-500

# mu_r = [0.1, 1, 10, 100]  # Regularization parameter
# mu_r = [1, 10, 50, 75, 100]  # Regularization parameter

########################################################################################################################
rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
R_V = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
rsv = 0.5
bck = 0.1  # Background water saturation time-lapse change
x_range = 157  # This is the # of grid cells in the x-direction
y_range = 159  # This is the # of grid cells in the y-direction
z_range = 5  # This is the # of grid cells in the z-direction
N = x_range_data * y_range_data  # This is the # number of data points/measurement stations
M = x_range * y_range * z_range  # This is the # number of model parameters


########################################################################################################################
########################################################################################################################
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


def data_weighting_matrix(flr, percent, observed_data):
    v = np.zeros((N, N))
    observed_data = observed_data.flatten()
    for j in range(N):
        for i in range(N):
            if i == j:
                v[j, i] = 1 / (flr * 1e-8 + (percent * observed_data[i]) * 1e-8)
    return v


def depth_weighting_matrix(dimension, z, z0, beta):
    w_model = np.zeros((dimension, dimension))
    for i in range(dimension):
        for j in range(dimension):
            if i == j:
                w_model[j, i] = pow((z[i] + z0), -beta / 2)
                # w_model[j, i] = pow((z[i] + z0), 1)
    return w_model


def grav_data_plot(dt, iteration, run_name, reg):
    X = np.linspace(0, 13213.48, 157)
    Y = np.linspace(0, 14173.40, 159)
    plt.style.use('classic')
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')
    Z = gaussian_filter(dt, sigma=0, mode='reflect')
    plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
               cmap=plt.cm.get_cmap('jet', 1000),
               # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
               interpolation='nearest', origin='upper')
    clb = plt.colorbar()
    tick_locator = ticker.MaxNLocator(nbins=7)
    clb.locator = tick_locator
    clb.update_ticks()
    clb.ax.set_title('$(\u03BC$Gal)', fontsize=15)
    ax.set_ylabel("Y(m)", labelpad=15)
    plt.title("" + str(float(delta_SW_starting_model)) + "_" + str(iteration) + "_" + str(reg) + "_" + str(run_name))
    label0 = str(delta_SW_starting_model).replace(".", ",")
    label1 = str(reg).replace(".", ",")
    label = "gv_result_" + label0 + "_" + str(iteration) + "_" + label1 + "_" + str(run_name)
    plt.xlabel("X [$km$]")
    plt.ylabel("Y [$km$]")
    plt.clim(0, 13)
    ax = plt.gca()
    ax.patch.set_visible(False)
    fig.patch.set_visible(False)
    plt.show()
    fig.savefig(label, transparent=True)
    plt.close(fig)


########################################################################################################################
#         Conjugate Gradient (CG) Inversion Algorithm

def conjugate_gradient_inversion(dsw0, mu, g, v, w, d_observed, d_initial):
    information = "None"
    d_initial = d_initial.flatten() * 1e-8
    d_observed = d_observed.flatten() * 1e-8
    misfit_r = (1 / 2) * norm(v * (np.subtract(d_initial, d_observed)))
    # model_norm_initial = (1 / 2) * norm(w * (dsw0 - ref_model))
    model_norm_initial = (1 / 2) * norm(w * dsw0)
    initial_total_misfit = misfit_r + model_norm_initial
    old_misfit = 0
    g_transpose = g.transpose()  # Please look at slide 360
    v_transpose = v.transpose()  # Please look at slide 360
    # w_transpose = w.transpose()  # Please look at slide 360

    # A = (((g_transpose @ v_transpose) @ v) @ g) + (mu * (w_transpose @ w))
    A = (((g_transpose @ v_transpose) @ v) @ g) + (mu * (pow(w, 2)))
    b = (((g_transpose @ v_transpose) @ v) @ d_observed)
    # b = (((g_transpose @ v_transpose) @ v) @ d_observed) + (mu * (pow(w, 2)) @ ref_model)
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
        if 1600 < misfit_r < 3100:  # Target misfit tolerance
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

        # modeled_data = grav_mod_coarse_real(dsw_k, count)
        modeled_data = gauss_fft_fwm(dsw_k, count)
        ################################################################################################################
        np.savetxt('recovered_model_temp_' + str(fig_num) + "_" + str(count) + "_" + str(mu) + '.txt', dsw_k,
                   delimiter=' ', header='')
        np.savetxt('recovered_data_temp_' + str(fig_num) + "_" + str(count) + "_" + str(mu) + '.txt', modeled_data,
                   delimiter=' ', header='')
        grav_data_plot(np.real(modeled_data), count, fig_num, mu)
        ################################################################################################################
        modeled_data = modeled_data.flatten() * 1e-8
        misfit_temp = (1 / 2) * norm(v * (np.subtract(modeled_data, d_observed)))
        model_norm_temp = (1 / 2) * norm(w * (dsw_k))
        total_misfit = misfit_temp + model_norm_temp
        misfit_r = misfit_temp

        if old_misfit == misfit_r:
            information = "The old and new misfits are identical!"
            del A, b
            break

        if np.abs(old_misfit - misfit_r) < 1:
            information = "np.abs(old_misfit - misfit_r) < 1 is achieved"
            del A, b
            break

        old_misfit = misfit_r
        count += 1
        iteration.append(count)
        misfit_iteration.append(total_misfit)

        print("This is the value of the data misfit for run " + str(fig_num) + " : ", misfit_r)
        print("This is the value of the model norm for run " + str(fig_num) + " : ", model_norm_temp)
        print("This is the value of the total misfit for run " + str(fig_num) + " : ", total_misfit)
        ################################################################################################################
        # fig1 = plt.figure(figsize=(6.4, 4.8))
        # ax = fig1.add_subplot(111)
        # ax.set_aspect('auto')
        # plt.plot(iteration, misfit_iteration)
        # plt.xlabel("Number of Iterations")
        # plt.ylabel("Total Misfit ($\phi_m + \phi_d$)")
        # plt.title("inversion_iteration_" + str(fig_num))
        # ax = plt.gca()
        # fig1.savefig("inversion_iteration_" + str(fig_num))
        ################################################################################################################
    plt.plot(iteration, misfit_iteration, label="\u03bc= " + str(mu))
    plt.xlabel("Number of Iterations")
    plt.ylabel("Total Misfit ($\phi_m + \phi_d$)")
    plt.yscale("log")
    plt.title("inversion_iteration_" + str(fig_num) + "_Z0=" + str(z0))
    ax = plt.gca()
    legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    plt.legend(loc="upper right")
    # Put a nicer background color on the legend.
    legend.get_frame()
    plt.savefig("inversion_iteration_t_" + str(fig_num))
    ####################################################################################################################
    # model_norm = (1 / 2) * norm(w * (dsw_k - ref_model))
    model_norm = (1 / 2) * norm(w * dsw_k)
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
    # dswat_start_model = dws_start_model(delta_SW_starting_model)
    dswat_start_model = dws_start_model(delta_SW_starting_model)
    np.savetxt('start_model_' + str(fig_num) + '.txt', dswat_start_model, delimiter=' ',
               header=' ')
    print("This is the run number:", fig_num)
    print("This is the initial guess in water saturation:", delta_SW_starting_model)
    print("This is the regularization parameter:", mu_r)
    print("This is z0:", z0)
    print("This is the data grid:", x_range_data, y_range_data)
    print("This is the time-step :", timestep)
    print("The FWM of the initial start model is running...")
    # data_0 = grav_mod_coarse_initial(dswat_start_model)
    # np.savetxt('grav_data_' + str(timestep) + '_start_model_Wisting_' + str(fig_num) + '_cluster.txt', data_0, delimiter=' ',
    #            header=' ')
    print("The FWM of the initial start model is done!")

    # ####################################################################################################################
    # #         OBSERVED MODEL
    # ####################################################################################################################
    if timestep == 32:
        dswat_observed = np.subtract(sat32, sat27)

    else:
        dswat_observed = np.subtract(sat57, sat27)

    if dt_grid_resolution == "c":
        data_0 = np.loadtxt('grav_data_start_model_dt_grid.txt', delimiter=' ')
    else:
        data_0 = np.loadtxt('grav_data_32_start_model_Wisting_129_cluster.txt', delimiter=' ')

    # ####################################################################################################################
    # # Adding Gaussian noise to the synthetic observed data
    # ####################################################################################################################
    mud = 0
    sigma = 0.02
    max_value = 50
    min_value = 0
    noise = np.random.normal(mud, sigma, data_observed.shape)
    noise_data_observed = data_observed + noise
    # noise_data_observed = data_observed
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
    # np.savetxt('grav_data_coarse_grid_coarse_dtpoint_grid_noisy.txt', noise_data_observed, delimiter=' ', header='Noisy')
    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################
    # Conjugate Gradient Inversion Analysis
    print("Conjugate Gradient Inversion Analysis starts now!")
    ####################################################################################################################
    flor = 0.001
    percent = 0.02
    beta = 2
    z = [600, 700, 800]
    # mu_r = 4.6
    # mu_r = [0.1, 1, 10, 100, 1000, 10000, 100000]
    # mu_r = np.linspace(3000, 5000, NUMBER_OF_ITERATIONS)
    misfit_d = np.zeros(NUMBER_OF_ITERATIONS)
    model_norm_d = np.zeros(NUMBER_OF_ITERATIONS)
    V_r = data_weighting_matrix(flor, percent, noise_data_observed)
    W_r = depth_weighting_matrix(M, dpt, z0, beta)

    # Iterating through different and multiple values of the regularization parameter!

    for i in range(NUMBER_OF_ITERATIONS):
        print("#######################################################################################################")
        print("Inversion iteration from scratch: ", i, " is running, for mu: ", mu_r[i])

        recovered_model_sol, mu_sol, misfit_sol, model_norm_sol, iteration_number_sol, total_misfit_sol, solution_model = \
            None, None, None, None, None, None, None

        recovered_model, initial_model, misfit, model_norm_r, direction, iteration_number, misfit_iteration_result = \
            conjugate_gradient_inversion(dswat_start_model, mu_r[i], gv_kernel, V_r, W_r, noise_data_observed, data_0)

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

        if 100 < misfit_d[i] < 4000:
            solution_model = recovered_model
            mu_sol = mu_r[i]
            misfit_sol = misfit
            model_norm_sol = model_norm_r
            iteration_number_sol = iteration_number_d
            total_misfit_sol = total_misfit_d
    #
    # ####################################################################################################################
    # Plotting the data misfit as a function of the regularization parameter
    plt.plot(mu_r, misfit_d)
    plt.plot(mu_sol, misfit_sol, "*", color="tab:blue", markersize=7)
    plt.axhline(y=9, color='orange', linestyle='-')
    plt.xlabel("Regularization parameter " u"(\u03bc)")
    plt.ylabel("Data misfit ($\phi_d$)")
    # plt.ylim([0, 50])
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("data_misfit_" + str(fig_num))
    # ####################################################################################################################
    # Plotting the model norm as a function of the regularization parameter
    plt.plot(mu_r, model_norm_d)
    plt.plot(mu_sol, model_norm_sol, "*", color="tab:blue", markersize=7)
    plt.xlabel("Regularization parameter " u"(\u03bc)")
    plt.ylabel("Model norm ($\phi_m$)")
    # plt.ylim([0, 50])
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("model_norm_" + str(fig_num))
    # ####################################################################################################################
    # Plotting the total misfit as a function of the iteration number
    plt.plot(mu_r, total_misfit_sol)
    plt.xlabel("Regularization parameter " u"(\u03bc)")
    plt.ylabel("Total Misfit ($\phi_m + \phi_d$)")
    plt.savefig("total_misfit_" + str(fig_num))
    # ####################################################################################################################
    # ####################################################################################################################
    print(
        "*************************************************************************************************************")

    print("These are the solutions:", mu_sol, misfit_sol, model_norm_sol)
    print(
        "*************************************************************************************************************")
    # ####################################################################################################################
    np.savetxt('solution_details' + str(fig_num) + '.txt', (mu_sol, misfit_sol, model_norm_sol), delimiter=' ',
               header='These are the solutions: mu_sol, data_misfit_sol, model_norm_sol ')
    print("Forward Modeling computation of the recovered model...")
    # recovered_data = grav_mod_coarse_final(solution_model)
    txt = "Last step in the inversion is running ..."
    recovered_data = gauss_fft_fwm(solution_model, txt)

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
# kernel = gravity_kernel()
# alpha = np.linspace(0, 13213.48, 157)
# beta = np.linspace(0, 14173.40, 159)
# z = dpt
# A = kernel[100][:]
#
# x = np.linspace(0, 13213.48, 55)
# y = np.linspace(0, 14173.40, 57)
# h = 400
# w = np.zeros(5)
# kernel = np.zeros(M)
# z0 = -400
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
#

# w = w ** 2
# plt.plot(w, z, label="Depth weighting")
# plt.plot(kernel[100][:], z, label="Kernel decay with depth")
# ax = plt.gca()
# ax.xaxis.tick_top()
# plt.gca().invert_yaxis()
# legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
# legend.get_frame()
# plt.xscale(value="log")
# plt.show()
#
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
