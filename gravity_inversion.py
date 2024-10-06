import matplotlib.pyplot as plt
import numpy as np
from numpy import *
# import cv2
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from Grav_Modeling_coarse_grid import grav_mod_coarse
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
delta_SW_starting_model = 0.1  # too optimistic guess
int_min = 0.25
int_max = 0.5
rsv = 0.5
bck = 0.1  # Background water saturation time-lapse change
# delta_SW_observed_model = np.array([bck, bck, bck, bck, bck, bck, bck, bck, bck, rn.uniform(int_min, int_max),
#                                     rn.uniform(int_min, int_max), rn.uniform(int_min, int_max),
#                                     rn.uniform(int_min, int_max), rn.uniform(int_min, int_max),
#                                     rn.uniform(int_min, int_max), rn.uniform(int_min, int_max),
#                                     rn.uniform(int_min, int_max), rn.uniform(int_min, int_max), bck, bck, bck, bck, bck,
#                                     bck, bck, bck, bck])  # realistic guess of the observed density model

delta_SW_observed_model = np.array([bck, bck, bck, bck, bck, bck, bck, bck, bck, rsv, rsv, rsv, rsv, rsv, rsv, rsv,
                                    rsv, rsv, bck, bck, bck, bck, bck, bck, bck, bck,
                                    bck])  # realistic guess of the observed density model
phi_starting_model = np.ones(27) * 0.3  # porosity of the start model
cell_dpt_starting_model = np.array([600, 600, 600, 600, 600, 600, 600, 600, 600, 700, 700, 700, 700, 700, 700, 700, 700,
                                    700, 800, 800, 800, 800, 800, 800, 800, 800, 800])
dx_accum_starting_model = np.array([380, 450, 520, 380, 450, 520, 380, 450, 520, 380, 450, 520, 380, 450, 520, 380, 450,
                                    520, 380, 450, 520, 380, 450, 520, 380, 450, 520])
dy_accum_starting_model = np.array([380, 450, 520, 380, 450, 520, 380, 450, 520, 380, 450, 520, 380, 450, 520, 380, 450,
                                    520, 380, 450, 520, 380, 450, 520, 380, 450, 520])
dx_starting_model = np.ones(27) * 70  # dx of the starting model
dy_starting_model = np.ones(27) * 70  # dy of the starting model
dz_starting_model = np.ones(27) * 2  # dz of the starting model
x_range = 3  # This is the # of grid cells in the x-direction
y_range = 3  # This is the # of grid cells in the y-direction
z_range = 3  # This is the # of grid cells in the z-direction
x_max = 900
y_max = 900
z_max = -1200
x_scale_starting_model = 300
y_scale_starting_model = 300
N = 9  # This is the # number of data points/measurement stations
D = 3  # This is the # number of data points/measurement stations along each horizontal direction


########################################################################################################################
# This is a script that computes the starting rho0 density or model vector to be used in the inversion scheme
def density_start_model(dsw):
    rho0 = (np.ones(x_range * y_range * z_range)) * (
            (rho_water - rho_oil) * dsw)
    return rho0


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
                v[j, i] = 1 / (0.001 + (0.03 * observed_data[i, j]))
    misfit = (1 / 2) * norm(v * (modeled_data - observed_data))
    return misfit


def gravity_kernel(cell_dpt_input, dx_accum_input, dy_accum_input,
                   x_range_input, y_range_input, z_range_input, x_scale_input, y_scale_input):
    G_kernel = np.zeros((x_range_input * y_range_input, x_range_input * y_range_input * z_range_input))
    a = (6.67 * 1e-11) * (cell_dpt_input - 400)
    n = 0

    for y_step in range(0, y_range_input):
        for x_step in range(0, x_range_input):
            m = 0
            for k in range(z_range_input):
                for j in range(y_range_input):
                    for i in range(x_range_input):
                        r = np.sqrt((((x_step * x_scale_input) + (x_scale_input / 2)) - dx_accum_input[m]) ** 2 + (
                                ((y_step * y_scale_input) + (y_scale_input / 2)) - dy_accum_input[m]) ** 2)
                        d = pow((pow(r, 2) + pow(cell_dpt_input[m] - 400, 2)), 1.5)
                        G_kernel[n, m] = ((a[m]) / d)
                        m += 1
            n += 1
    return G_kernel


def data_weighting_matrix(dimension, floor, percent, observed_data):
    v = np.zeros((dimension, dimension))
    observed_data = observed_data.flatten()
    for i in range(dimension):
        for j in range(dimension):
            if i == j:
                v[j, i] = 1 / (0.001 * 1e-8 + (0.03 * observed_data[i]) * 1e-8)
                v[j, i] = 1 / (floor * 1e-8 + (percent * observed_data[i]) * 1e-8)
    return v


def depth_weighting_matrix(dimension, z, z0, beta):
    w = np.zeros((dimension, dimension))
    for i in range(dimension):
        for j in range(dimension):
            if i == j:
                w[j, i] = pow((z[i] + z0), -beta / 2)
    return w


########################################################################################################################
#         Conjugate Gradient (CG) Inversion Algorithm

def conjugate_gradient_inversion(rho0, mu, g, v, w, d_observed, d_initial):
    d_initial = d_initial.flatten() * 1e-8
    d_observed = d_observed.flatten() * 1e-8
    misfit_r = (1 / 2) * norm(v * (np.subtract(d_initial, d_observed)))
    model_norm_initial = (1 / 2) * norm(w * rho0)
    initial_total_misfit = misfit_r + model_norm_initial
    old_misfit = 0
    g_transpose = g.transpose()  # Please look at slide 360
    v_transpose = v.transpose()  # Please look at slide 360
    w_transpose = w.transpose()  # Please look at slide 360

    A = (((g_transpose @ v_transpose) @ v) @ g) + (mu * (w_transpose @ w))
    b = (((g_transpose @ v_transpose) @ v) @ d_observed)
    r0 = (A @ rho0) - b
    p0 = -r0
    rho_k = rho0
    rk = r0
    pk = p0

    count = 0
    iteration = [0]
    misfit_iteration = [initial_total_misfit]

    while rk.all() != 0:
        if 3 < misfit_r < 9:  # Target misfit tolerance
            break

        alpha = (rk.transpose() @ rk) / ((pk.transpose() @ A) @ pk)
        rho_new = rho_k + (alpha * pk)
        r_new = rk + (alpha * (A @ pk))
        beta_r = (r_new.transpose() @ r_new) / (rk.transpose() @ rk)
        p_new = -r_new + (beta_r * pk)

        rho_k = rho_new
        rk = r_new
        pk = p_new

        modeled_data = grav_mod_coarse(rho_k, cell_dpt_starting_model, phi_starting_model, dx_starting_model,
                                       dy_starting_model, dz_starting_model, dx_accum_starting_model,
                                       dy_accum_starting_model, x_range, y_range, z_range, x_scale_starting_model,
                                       y_scale_starting_model)
        modeled_data = modeled_data.flatten() * 1e-8
        misfit_temp = (1 / 2) * norm(v * (np.subtract(modeled_data, d_observed)))
        model_norm_temp = (1 / 2) * norm(w * (rho_k))
        total_misfit = misfit_temp + model_norm_temp
        misfit_r = misfit_temp

        if old_misfit == misfit_r:
            break
        old_misfit = misfit_r
        count += 1
        iteration.append(count)
        misfit_iteration.append(total_misfit)

        print("Still iterating!")
        print("This is the value of the misfit:", misfit_r)

    model_norm = (1 / 2) * norm(w * (rho_k))
    print("The Conjugate Gradient Inversion is successfully done Alhamdulilah")
    return rho_k, rho0, misfit_r, model_norm, pk, iteration, misfit_iteration


def gravity_inversion_execution():
    ####################################################################################################################

    #         FUNCTION EXECUTIONS

    ####################################################################################################################
    i = 0
    # Starting model for the inversion
    ####################################################################################################################
    #         STARTING MODEL
    ####################################################################################################################
    print(
        "WELCOME TO THE GRAVITY INVERSION SCRIPT, FROM HERE YOU WILL BE PROVIDED WITH rho_0, rho_observed, data_0 AND "
        "data_observed ")
    ####################################################################################################################
    rho_starting_model = density_start_model(delta_SW_starting_model)
    np.savetxt('rho_starting_model.txt', rho_starting_model, delimiter=' ', header=' ')
    print("FWM for iteration: ", i, " is running...")
    data_0 = grav_mod_coarse(rho_starting_model, cell_dpt_starting_model, phi_starting_model, dx_starting_model,
                             dy_starting_model,
                             dz_starting_model, dx_accum_starting_model, dy_accum_starting_model, x_range, y_range,
                             z_range,
                             x_scale_starting_model, y_scale_starting_model)
    np.savetxt('grav_inv_dummy_starting_model.txt', data_0, delimiter=' ', header=' ')
    ####################################################################################################################
    #         OBSERVED MODEL
    ####################################################################################################################
    rho_observed = density_start_model(delta_SW_observed_model)
    np.savetxt('rho_observed.txt', rho_observed, delimiter=' ', header=' ')
    data_observed = grav_mod_coarse(rho_observed, cell_dpt_starting_model, phi_starting_model, dx_starting_model,
                                    dy_starting_model, dz_starting_model, dx_accum_starting_model,
                                    dy_accum_starting_model,
                                    x_range, y_range, z_range, x_scale_starting_model, y_scale_starting_model)
    np.savetxt('grav_inv_dummy_observed_model.txt', data_observed, delimiter=' ', header=' ')
    ####################################################################################################################
    # Adding Gaussian noise to the synthetic observed data
    ####################################################################################################################
    mud = 0
    sigma = 0.005
    noise = np.random.normal(mud, sigma, data_observed.shape)
    noise_data_observed = data_observed + noise
    data_number = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    d1 = data_observed.flatten()
    d2 = noise_data_observed.flatten()
    # plt.plot(data_number, d1, "+", color='black', label="Accurate data")
    # plt.plot(data_number, d2, "o", fillstyle='none', color='black', label="Noisy data")
    # plt.xlabel("data number")
    # plt.ylabel("Gravity data $(\u03BC$Gal)")
    # ax = plt.gca()
    # legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    # plt.legend(loc="upper right")
    # Put a nicer background color on the legend.
    # legend.get_frame()
    # plt.show()
    ####################################################################################################################
    # Conjugate Gradient Inversion Analysis
    ####################################################################################################################
    floor = 0.001
    percent = 0.03
    beta = 2
    # z0 = -400
    z0 = -375
    z = [600, 700, 800]
    # mu_r = 4.6
    NUMBER_OF_ITERATIONS = 1000
    mu_r = np.linspace(0.0001, 100, NUMBER_OF_ITERATIONS)
    misfit_d = np.zeros(NUMBER_OF_ITERATIONS)
    model_norm_d = np.zeros(NUMBER_OF_ITERATIONS)
    G_r = gravity_kernel(cell_dpt_starting_model, dx_accum_starting_model, dy_accum_starting_model, x_range, y_range,
                         z_range, x_scale_starting_model, y_scale_starting_model)
    V_r = data_weighting_matrix(N, floor, percent, noise_data_observed)
    W_r = depth_weighting_matrix(N * D, cell_dpt_starting_model, z0, beta)

    # Iterating through different and multiple values of the regularization parameter!

    recovered_model_sol, mu_sol, misfit_sol, model_norm_sol, iteration_number_sol, total_misfit_sol = \
        None, None, None, None, None, None

    for i in range(NUMBER_OF_ITERATIONS):
        recovered_model, initial_model, misfit, model_norm_r, direction, iteration_number, misfit_iteration_result = \
            conjugate_gradient_inversion(rho_starting_model, mu_r[i], G_r, V_r, W_r, noise_data_observed, data_0)
        print("This is the initial model", rho_starting_model)
        print("This is the observed model", rho_observed)
        print("This is the recovered model", recovered_model)
        print("This is the data misft:", misfit)
        print("This is the model norm:", model_norm_r)
        print("This is the regularization parameter:", mu_r[i])
        misfit_d[i] = misfit
        model_norm_d[i] = model_norm_r
        total_misfit_d = misfit_iteration_result
        iteration_number_d = iteration_number

        if 8.5 < misfit_d[i] < 10:
            recovered_model_sol = recovered_model
            mu_sol = mu_r[i]
            misfit_sol = misfit
            model_norm_sol = model_norm_r
            iteration_number_sol = iteration_number_d
            total_misfit_sol = total_misfit_d

    ####################################################################################################################
    # Plotting the data misfit as a function of the regularization parameter
    fig0 = plt.figure(figsize=(6.4, 4.8))
    plt.plot(mu_r, misfit_d)
    plt.plot(mu_sol, misfit_sol, "*", color="tab:blue", markersize=7)
    plt.axhline(y=9, color='orange', linestyle='-')
    plt.xlabel("Regularization parameter " u"(\u03bc)")
    plt.ylabel("Data misfit ($\phi_d$)")
    # plt.ylim([0, 50])
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    fig0.savefig("Data_misfit_vs_mu.png")
    ####################################################################################################################
    # Plotting the model norm as a function of the regularization parameter
    fig1 = plt.figure(figsize=(6.4, 4.8))
    plt.plot(mu_r, model_norm_d)
    plt.plot(mu_sol, model_norm_sol, "*", color="tab:blue", markersize=7)
    plt.xlabel("Regularization parameter " u"(\u03bc)")
    plt.ylabel("Model norm ($\phi_m$)")
    # plt.ylim([0, 50])
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    fig1.savefig("Model_norm_vs_mu.png")
    ####################################################################################################################
    # Plotting the total misfit as a function of the iteration number
    fig2 = plt.figure(figsize=(6.4, 4.8))
    plt.plot(iteration_number_sol, total_misfit_sol)
    # plt.plot(mu_sol, model_norm_sol, "*", color="tab:blue", markersize=7)
    plt.xlabel("Number of Iterations")
    plt.ylabel("Total Misfit ($\phi_m + \phi_d$)")
    # plt.show()
    fig2.savefig("Total_misfit_vs_nr_of_iterations.png")
    ####################################################################################################################
    ####################################################################################################################
    print(
        "*************************************************************************************************************")

    print("These are the solutions:", mu_sol, misfit_sol, model_norm_sol)
    print(
        "*************************************************************************************************************")
    ####################################################################################################################
    recovered_data = grav_mod_coarse(recovered_model_sol, cell_dpt_starting_model, phi_starting_model,
                                     dx_starting_model, dy_starting_model,
                                     dz_starting_model, dx_accum_starting_model, dy_accum_starting_model, x_range,
                                     y_range, z_range,
                                     x_scale_starting_model, y_scale_starting_model)

    return rho_starting_model, rho_observed, recovered_model_sol, recovered_data, data_observed

########################################################################################################################
########################################################################################################################
#         Data Misfit
########################################################################################################################
########################################################################################################################
# Depth weighting analysis
########################################################################################################################
beta = 2
z0 = -375
z = [600, 700, 800]
w = depth_weighting_matrix(N * D, cell_dpt_starting_model, z0, beta)
w_transpose = w.transpose()
A = (10 * (w_transpose @ w))
alpha = [380, 450, 520]
beta = [380, 450, 520]
z = [600, 700, 800]
h = 400
w = np.zeros(3)
kernel = np.zeros(3)
z0 = -375
for i in range(3):
    w[i] = (1 / (z[i] + z0))
    r = np.sqrt((450 - alpha[i]) ** 2 + (
            450 - beta[i]) ** 2)
    d = pow((pow(r, 2) + pow(z[i] - 400, 2)), 1.5)
    kernel[i] = (z[i] - h) / d
    # kernel[i] = 1 / ((z[i] - h) ** 2)

w = w ** 2
plt.plot(w, z, label="Depth weighting")
plt.plot(kernel, z, label="Kernel decay with depth")
ax = plt.gca()
ax.xaxis.tick_top()
plt.gca().invert_yaxis()
legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
legend.get_frame()
plt.show()

plt.plot(w, z, label="Depth weighting")
plt.plot(kernel, z, label="Kernel decay with depthh")
ax = plt.gca()
plt.xscale(value="log")
ax.xaxis.tick_top()
plt.gca().invert_yaxis()
legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
plt.xlabel("Kernel value")
plt.ylabel("Depth (m)")
legend.get_frame()
plt.show()
