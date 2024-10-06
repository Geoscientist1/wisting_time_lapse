from numpy import *
import numpy as np
import time

start_time = time.time()

# from files_readme import celltopdpt, phi, cell_vol
# from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
# from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
# from files_readme import phi, cel_cen_dpt, dx, dy, dz, dx_accum, dy_accum, A, B

# dx = np.loadtxt('dx_coarse.txt', delimiter=' ')
# dy = np.loadtxt('dy_coarse.txt', delimiter=' ')
# dz = np.loadtxt('dz_coarse.txt', delimiter=' ')
# dv = np.loadtxt('dv_coarse.txt', delimiter=' ')
# dx_accum = np.loadtxt('dx_accum_coarse.txt', delimiter=' ')
# dy_accum = np.loadtxt('dy_accum_coarse.txt', delimiter=' ')
# phi = np.loadtxt('phi_coarse.txt', delimiter=' ')
cel_cen_dpt = np.loadtxt('cel_cen_dpt_coarse.txt', delimiter=' ')
dpt_cell = np.loadtxt('dpt_cell_initial_model_inversion.txt', delimiter=' ')
dx_cell = np.loadtxt('dx_cell_initial_model_inversion.txt', delimiter=' ')
dy_cell = np.loadtxt('dy_cell_initial_model_inversion.txt', delimiter=' ')
dz_cell = np.loadtxt('dz_cell_initial_model_inversion.txt', delimiter=' ')

# rho27_27 = np.loadtxt('rho27_27_coarse.txt', delimiter=' ')
# rho27_29 = np.loadtxt('rho27_29_coarse.txt', delimiter=' ')
# rho27_31 = np.loadtxt('rho27_31_coarse.txt', delimiter=' ')
# rho27_33 = np.loadtxt('rho27_33_coarse.txt', delimiter=' ')
# rho27_35 = np.loadtxt('rho27_35_coarse.txt', delimiter=' ')
# rho27_39 = np.loadtxt('rho27_39_coarse.txt', delimiter=' ')
# rho27_45 = np.loadtxt('rho27_45_coarse.txt', delimiter=' ')
# rho27_49 = np.loadtxt('rho27_49_coarse.txt', delimiter=' ')
# rho27_53 = np.loadtxt('rho27_53_coarse.txt', delimiter=' ')
# rho27_57 = np.loadtxt('rho27_57_coarse.txt', delimiter=' ')

sat27 = np.loadtxt('sat27_coarse.txt', delimiter=' ')
sat32 = np.loadtxt('sat32_coarse.txt', delimiter=' ')
sat42 = np.loadtxt('sat42_coarse.txt', delimiter=' ')
sat57 = np.loadtxt('sat57_coarse.txt', delimiter=' ')
# ref_model = np.loadtxt('reference_model.txt', delimiter=' ')
# sat31 = np.loadtxt('sat31_coarse.txt', delimiter=' ')
# sat33 = np.loadtxt('sat33_coarse.txt', delimiter=' ')
# sat35 = np.loadtxt('sat35_coarse.txt', delimiter=' ')
# sat39 = np.loadtxt('sat39_coarse.txt', delimiter=' ')
# sat45 = np.loadtxt('sat45_coarse.txt', delimiter=' ')
# sat49 = np.loadtxt('sat49_coarse.txt', delimiter=' ')
# sat53 = np.loadtxt('sat53_coarse.txt', delimiter=' ')
# sat57 = np.loadtxt('sat57_coarse.txt', delimiter=' ')

# soil_27 = np.loadtxt('soil_27_coarse.txt', delimiter=' ')
# soil_57 = np.loadtxt('soil_57_coarse.txt', delimiter=' ')
# soil_31 = np.loadtxt('soil_31_coarse.txt', delimiter=' ')
# soil_33 = np.loadtxt('soil_33_coarse.txt', delimiter=' ')
# soil_35 = np.loadtxt('soil_35_coarse.txt', delimiter=' ')
# soil_39 = np.loadtxt('soil_39_coarse.txt', delimiter=' ')
# soil_45 = np.loadtxt('soil_45_coarse.txt', delimiter=' ')
# soil_49 = np.loadtxt('soil_49_coarse.txt', delimiter=' ')
# soil_53 = np.loadtxt('soil_53_coarse.txt', delimiter=' ')
# soil_57 = np.loadtxt('soil_57_coarse.txt', delimiter=' ')

rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
########################################################################################################################
########################################################################################################################
x_range = 157  # This is the # of grid cells in the x-direction
y_range = 159  # This is the # of grid cells in the y-direction
z_range = 5  # This is the # of grid cells in the z-direction
M = x_range * y_range * z_range  # This is the # number of model parameters
k_range = 101
k_range_new = 5
########################################################################################################################
########################################################################################################################
dt_grid_resolution = "c"
timestep = 57
if dt_grid_resolution == "f":
    x_range_data = 157  # This is the # of grid data points  in the x-direction
    y_range_data = 159  # This is the # of grid data points  in the y-direction
    XSCALE = 84
    YSCALE = 89
    gv_kernel = np.loadtxt('gravity_kernel_157_159.txt', delimiter=' ')
    if timestep == 32:
        data_observed = np.loadtxt('grav_data_32_coarse.txt', delimiter=' ')
    else:
        data_observed = np.loadtxt('grav_data_57_coarse.txt', delimiter=' ')
elif dt_grid_resolution == "c":
    x_range_data = 55  # This is the # of grid data points  in the x-direction
    y_range_data = 57  # This is the # of grid data points  in the y-direction
    XSCALE = 240
    YSCALE = 249
    gv_kernel = np.loadtxt('gravity_kernel_55_57.txt', delimiter=' ')
    # gv_kernel = np.loadtxt('gravity_kernel_55_57_revised.txt', delimiter=' ')
    if timestep == 32:
        data_observed = np.loadtxt('grav_data_32_coarse_dt_grid.txt', delimiter=' ')
    else:
        data_observed = np.loadtxt('grav_data_57_coarse_dt_grid.txt', delimiter=' ')


########################################################################################################################
########################################################################################################################


# Now we account for the different water saturations corresponding to the different times from year 2027 all the way
# to year 2039
########################################################################################################################
########################################################################################################################
# This is the FWM engine used with the Wisting model
########################################################################################################################
# def grav_mod_coarse(rho_input, cell_dpt_input, phi_input, dx_input, dy_input, dz_input, dx_accum_input, dy_accum_input,
#                     x_range_input, y_range_input, z_range_input, x_scale_input, y_scale_input):


########################################################################################################################
# THIS IS THE OLD FORWARD MODELING ENGINE ############################################################################
########################################################################################################################
########################################################################################################################
def grav_mod_coarse_real(dsw, f):
    delta_g = np.zeros((y_range_data, x_range_data))
    delta_rho = dsw * (rho_water - rho_oil)
    a = delta_rho * phi
    b = dx * dy * dz  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    n = 0
    print("Forward Modeling computation...")
    print("Inversion iteration: " + str(f) + " is running...")
    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            m = 0
            for k in range(k_range_new):
                for j in range(159):
                    for i in range(157):
                        delta_g[y_step, x_step] += ((c[m]) * gv_kernel[n, m])
                        m += 1
            n += 1
    return delta_g * 1e8


########################################################################################################################
########################################################################################################################
########################################################################################################################
def gravity_kernel():
    G_kernel = np.zeros((x_range_data * y_range_data, M))
    a = G * (cel_cen_dpt - 400)
    dx_accumulated = np.linspace(0, 13213.48, 157)
    dy_accumulated = np.linspace(0, 14173.40, 159)
    n = 0
    print("Kernel computation...")
    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            print("y", y_step)
            print("x", x_step)
            m = 0
            for k in range(z_range):
                for j in range(y_range):
                    for i in range(x_range):
                        r = np.sqrt(
                            ((x_step * XSCALE) - dx_accumulated[i]) ** 2 + ((y_step * YSCALE) - dy_accumulated[j]) ** 2)
                        d = pow((pow(r, 2) + pow(cel_cen_dpt[m] - 400, 2)), 1.5)
                        G_kernel[n, m] = ((a[m]) / d)
                        m += 1
            n += 1
    return G_kernel


def gravity_kernel_revised():
    G_kernel = np.zeros((x_range_data * y_range_data, M))
    a = G * (dpt_cell - 400)
    dx_accumulated = np.linspace(0, 13213.48, 157)
    dy_accumulated = np.linspace(0, 14173.40, 159)
    n = 0
    print("Revised kernel computation...")
    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            print("y", y_step)
            print("x", x_step)
            m = 0
            for k in range(z_range):
                for j in range(y_range):
                    for i in range(x_range):
                        r = np.sqrt(
                            ((x_step * XSCALE) - dx_accumulated[i]) ** 2 + ((y_step * YSCALE) - dy_accumulated[j]) ** 2)
                        d = pow((pow(r, 2) + pow(dpt_cell[m] - 400, 2)), 1.5)
                        G_kernel[n, m] = ((a[m]) / d)
                        m += 1
            n += 1
    return G_kernel


########################################################################################################################
########################################################################################################################
########################################################################################################################
def grav_mod_coarse_final(sat_1, sat_2):
    delta_g = np.zeros((y_range_data, x_range_data))
    rho_1 = (sat_1 * rho_water) + ((1 - sat_1) * rho_oil)
    rho_2 = (sat_2 * rho_water) + ((1 - sat_2) * rho_oil)
    delta_rho = rho_2 - rho_1
    n = 0
    a = delta_rho * phi
    b = dx * dy * dz  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            m = 0
            for k in range(k_range_new):
                for j in range(159):
                    for i in range(157):
                        # r = np.sqrt(
                        #     ((x_step * XSCALE) - dx_accumulated[i]) ** 2 + ((y_step * YSCALE) - dy_accumulated[j]) ** 2)
                        # d = pow((pow(r, 2) + pow(cel_cen_dpt[m] - 400, 2)), 1.5)
                        # delta_g[y_step, x_step] += ((c[m]) / d)
                        delta_g[y_step, x_step] += ((c[m]) * gv_kernel[n, m])
                        m += 1
            n += 1
    return delta_g * 1e8


########################################################################################################################
# THIS IS THE FORWARD MODELING ENGINE FOR INSIDE THE INVERSION #########################################################
########################################################################################################################
########################################################################################################################
def grav_mod_coarse_inversion(dsw, f):
    delta_g = np.zeros((y_range_data, x_range_data))

    delta_rho = dsw * (rho_water - rho_oil)
    phi_cell = 0.3
    dx_accumulated = np.linspace(0, 13213.48, 157)
    dy_accumulated = np.linspace(0, 14173.40, 159)
    a = G * (dpt_cell - 400) * delta_rho * phi_cell
    b = dx_cell * dy_cell * dz_cell  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    print("Forward Modeling computation...")
    print("Inversion iteration: " + str(f) + " is running...")
    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            # print("y", y_step)
            # print("x", x_step)
            m = 0
            for k in range(k_range_new):
                for j in range(159):
                    for i in range(157):
                        r = np.sqrt(
                            ((x_step * XSCALE) - dx_accumulated[i]) ** 2 + ((y_step * YSCALE) - dy_accumulated[j]) ** 2)
                        d = pow((pow(r, 2) + pow(dpt_cell[m] - 400, 2)), 1.5)
                        delta_g[y_step, x_step] += ((c[m]) / d)
                        m += 1
    return delta_g * 1e8


########################################################################################################################
# THIS IS THE FORWARD MODELING ENGINE FOR THE INITIAL INVERSION MODEL ##################################################
########################################################################################################################
########################################################################################################################
def grav_mod_coarse_initial(dsw):
    delta_g = np.zeros((y_range_data, x_range_data))
    delta_rho = dsw * (rho_water - rho_oil)
    phi_cell = 0.3
    dx_accumulated = np.linspace(0, 13213.48, 157)
    dy_accumulated = np.linspace(0, 14173.40, 159)
    a = G * (dpt_cell - 400) * delta_rho * phi_cell
    b = dx_cell * dy_cell * dz_cell  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            print("y", y_step)
            print("x", x_step)
            m = 0
            for k in range(k_range_new):
                for j in range(159):
                    for i in range(157):
                        r = np.sqrt(
                            ((x_step * XSCALE) - dx_accumulated[i]) ** 2 + ((y_step * YSCALE) - dy_accumulated[j]) ** 2)
                        d = pow((pow(r, 2) + pow(dpt_cell[m] - 400, 2)), 1.5)
                        delta_g[y_step, x_step] += ((c[m]) / d)
                        m += 1
    return delta_g * 1e8


########################################################################################################################
########################################################################################################################
########################################################################################################################
# THIS IS THE FORWARD MODELING ENGINE WHEN USING THE ###################################################################
# MAPPING MATRIX #######################################################################################################
########################################################################################################################
########################################################################################################################
# def grav_mod_coarse(sat_1, sat_2):
def grav_mod_coarse(dsw, val, porr, cell_depthh, dx_accumulatedd, dy_accumulatedd, xd, yd, dzz, f):
    delta_g = np.zeros((y_range_data, x_range_data))

    # rho_1 = (sat_1 * rho_water) + ((1 - sat_1) * rho_oil)
    # rho_2 = (sat_2 * rho_water) + ((1 - sat_2) * rho_oil)
    # delta_rho = rho_2 - rho_1
    delta_rho = dsw * (rho_water - rho_oil)
    a = G * (cell_depthh - 400) * delta_rho * porr
    b = xd * yd * dzz  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    dx_accumulated = np.linspace(0, 13213.48, 157)
    dy_accumulated = np.linspace(0, 14173.40, 159)
    print("Forward Modeling computation...")
    print("Inversion iteration: " + str(f) + " is running...")
    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            # print("Inversion iteration: " + str(f) + " is running!")
            # print("Forward Modeling computation")
            # print("y", y_step)
            # print("x", x_step)
            m = 0
            for k in range(val):
                r = np.sqrt(
                    ((x_step * XSCALE) - dx_accumulatedd[m]) ** 2 + ((y_step * YSCALE) - dy_accumulatedd[m]) ** 2)
                d = pow((pow(r, 2) + pow(cell_depthh[m] - 400, 2)), 1.5)
                delta_g[56 - y_step, x_step] += ((c[m]) / d)
                m += 1
    return delta_g * 1e8


########################################################################################################################
########################################################################################################################
########################################################################################################################
# THIS FUNCTION TO MAKE THE CELL DEPTH CORRECT OF THE INITIAL MODEL ####################################################
########################################################################################################################
########################################################################################################################
def initial_model_depthing():
    dpt = np.zeros(len(cel_cen_dpt))
    cell_dx = np.zeros(len(cel_cen_dpt))
    cell_dy = np.zeros(len(cel_cen_dpt))
    cell_dz = np.zeros(len(cel_cen_dpt))
    m = 0
    print("Initial model depthing is running")
    for k in range(5):
        vari = 650 + (k * 5)
        for j in range(159):
            for i in range(157):
                dpt[m] = vari
                cell_dx[m] = 74
                cell_dy[m] = 74
                cell_dz[m] = 3
                m += 1
    return dpt, cell_dx, cell_dy, cell_dz


########################################################################################################################
########################################################################################################################
########################################################################################################################
# THIS FUNCTION TO CONSTRUCT THE REFERENCE MODEL FOR THE INVERSION #####################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

def reference_model():
    for i in range(len(sat27)):
        if sat27[i] == 1.0:
            sat27[i] = 0.0
    return sat27


########################################################################################################################
# G_r = gravity_kernel_revised()
# np.savetxt('gravity_kernel_55_57_revised.txt', G_r, delimiter=' ', header='')
########################################################################################################################
# ref_model = reference_model()
# np.savetxt('reference_model.txt', ref_model, delimiter=' ', header='')
########################################################################################################################
# dpt, cell_dx, cell_dy, cell_dz = initial_model_depthing()
# np.savetxt('dpt_cell_initial_model_inversion.txt', dpt, delimiter=' ', header='')
# np.savetxt('dx_cell_initial_model_inversion.txt', cell_dx, delimiter=' ', header='')
# np.savetxt('dy_cell_initial_model_inversion.txt', cell_dy, delimiter=' ', header='')
# np.savetxt('dz_cell_initial_model_inversion.txt', cell_dz, delimiter=' ', header='')
########################################################################################################################
# gv27_27 = grav_mod_coarse(sat27, sat27, soil_27, soil_27)
# gv27_29 = grav_mod_coarse(sat27, sat29, soil_27, soil_29)
# gv27_31 = grav_mod_coarse(sat27, sat31, soil_27, soil_31)
# gv27_33 = grav_mod_coarse(sat27, sat33, soil_27, soil_33)
# gv27_35 = grav_mod_coarse(sat27, sat35, soil_27, soil_35)
# gv27_39 = grav_mod_coarse(sat27, sat39, soil_27, soil_39)
# gv27_45 = grav_mod_coarse(sat27, sat45, soil_27, soil_45)
# gv27_49 = grav_mod_coarse(sat27, sat49, soil_27, soil_49)
# gv27_53 = grav_mod_coarse(sat27, sat53, soil_27, soil_53)

# gv27_57 = grav_mod_coarse_real(sat27, sat57)
# np.savetxt('fast_FWM_fine_dt_grid.txt', gv27_57, delimiter=' ', header='Gravity data Coarse Grid for year 2057')

# gv = grav_mod_coarse_initial(0.5)
# np.savetxt('grav_data_start_model_dt_grid.txt', gv, delimiter=' ', header='Gravity data Coarse Grid for year 2057')


# end_time = time.time()
# execution_time = start_time - end_time
# print("Execution time:", execution_time)
