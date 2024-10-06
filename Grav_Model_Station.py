from numpy import *
import numpy as np
import time

start_time = time.time()

# from files_readme import celltopdpt, phi, cell_vol
from files_readme import sat27, sat29, sat31, sat32, sat33, sat35, sat37, sat39, sat41, sat42, sat43, sat45, sat47, \
    sat49, sat51, sat53, sat55, sat57
from files_readme import phi, cel_cen_dpt, dx, dy, dz, A, B

rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant

x_range_data = 157
y_range_data = 159
k_range = 101
k_range_new = 5

XSCALE = 84
YSCALE = 89

# Now we account for the different water saturations corresponding to the different times from year 2027 all the way
# to year 2039
################################################################################
################################################################################

# gv27_27 = np.loadtxt("grav_data_27.txt", delimiter=" ")


# The function calculating the gravity responses delta gz
def grav_mod(sat_1, sat_2):
    delta_g = np.zeros((y_range_data, x_range_data))
    fixed1, fixed2, fixed3, fixed4, fixed5 = None, None, None, None, None

    rho_1 = (sat_1 * rho_water) + ((1 - sat_1) * rho_oil)
    rho_2 = (sat_2 * rho_water) + ((1 - sat_2) * rho_oil)
    delta_rho = rho_2 - rho_1
    dx_accumulated = np.linspace(0, 13213.48, 157)
    dy_accumulated = np.linspace(0, 14173.40, 159)
    a = G * (cel_cen_dpt - 400) * delta_rho * phi
    b = dx * dy * dz  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b

    for y_step in range(0, y_range_data):
        for x_step in range(0, x_range_data):
            print("y", y_step)
            print("x", x_step)
            print("42")
            m = 0
            for k in range(k_range):
                for j in range(159):
                    for i in range(157):
                        r = np.sqrt(
                            ((x_step * XSCALE) - dx_accumulated[i]) ** 2 + ((y_step * YSCALE) - dy_accumulated[j]) ** 2)
                        d = pow((pow(r, 2) + pow(cel_cen_dpt[m] - 400, 2)), 1.5)
                        delta_g[y_step, x_step] += ((c[m]) / d)
                        m += 1

            if x_step * XSCALE == 3612 and y_step * YSCALE == 11303:
                fixed1 = delta_g[y_step, x_step]
            if x_step * XSCALE == 8400 and y_step * YSCALE == 9345:
                fixed2 = delta_g[y_step, x_step]
            if x_step * XSCALE == 9156 and y_step * YSCALE == 7209:
                fixed3 = delta_g[y_step, x_step]
            if x_step * XSCALE == 11256 and y_step * YSCALE == 9612:
                fixed4 = delta_g[y_step, x_step]
            if x_step * XSCALE == 12348 and y_step * YSCALE == 5518:
                fixed5 = delta_g[y_step, x_step]
    return delta_g * 1e8, fixed1 * 1e8, fixed2 * 1e8, fixed3 * 1e8, fixed4 * 1e8, fixed5 * 1e8


################################################################################
print('Nooooooooooooooooooooooooooooooooooooooooooooooooooow')


# gv27_27, fixed_1_27, fixed_2_27, fixed_3_27, fixed_4_27, fixed_5_27 = grav_mod(sat27, sat27)
# np.savetxt('grav_data_27_101_main.txt', gv27_27, delimiter=' ', header='Gravity data for year 2027')
# gv27_29, fixed_1_29, fixed_2_29, fixed_3_29, fixed_4_29, fixed_5_29 = grav_mod(sat27, sat29)
# np.savetxt('grav_data_29_101_main.txt', gv27_29, delimiter=' ', header='Gravity data for year 2029')
# gv27_31, fixed_1_31, fixed_2_31, fixed_3_31, fixed_4_31, fixed_5_31 = grav_mod(sat27, sat31)
# np.savetxt('grav_data_31_101_main.txt', gv27_31, delimiter=' ', header='Gravity data for year 2031')
gv27_42, fixed_1_42, fixed_2_42, fixed_3_42, fixed_4_42, fixed_5_42 = grav_mod(sat27, sat42)
np.savetxt('grav_data_42_101_main.txt', gv27_42, delimiter=' ', header='Gravity data for year 2042')
# gv27_35, fixed_1_35, fixed_2_35, fixed_3_35, fixed_4_35, fixed_5_35 = grav_mod(sat27, sat35)
# np.savetxt('grav_data_35_101_main.txt', gv27_35, delimiter=' ', header='Gravity data for year 2035')
# gv27_57, fixed_1_57, fixed_2_57, fixed_3_57, fixed_4_57, fixed_5_57 = grav_mod(sat27, sat57)
# np.savetxt('grav_data_57_101_main.txt', gv27_57, delimiter=' ', header='Gravity data for year 2057')
# gv27_39, fixed_1_39, fixed_2_39, fixed_3_39, fixed_4_39, fixed_5_39 = grav_mod(sat27, sat39)
# np.savetxt('grav_data_39_101_main.txt', gv27_39, delimiter=' ', header='Gravity data for year 2039')
# gv27_45, fixed_1_45, fixed_2_45, fixed_3_45, fixed_4_45, fixed_5_45 = grav_mod(sat27, sat45)
# np.savetxt('grav_data_45_101_main.txt', gv27_45, delimiter=' ', header='Gravity data for year 2045')
# gv27_49, fixed_1_49, fixed_2_49, fixed_3_49, fixed_4_49, fixed_5_49 = grav_mod(sat27, sat49)
# np.savetxt('grav_data_49_101_main.txt', gv27_49, delimiter=' ', header='Gravity data for year 2049')
# gv27_53, fixed_1_53, fixed_2_53, fixed_3_53, fixed_4_53, fixed_5_53 = grav_mod(sat27, sat53)
# np.savetxt('grav_data_53_101_main.txt', gv27_53, delimiter=' ', header='Gravity data for year 2053')
# gv27_57, fixed_1_57, fixed_2_57, fixed_3_57, fixed_4_57, fixed_5_57 = grav_mod(sat27, sat57)
# np.savetxt('grav_data_57_101_main.txt', gv27_57, delimiter=' ', header='Gravity data for year 2057')
########################################################################################################################
# fixed_27 = np.array([fixed_1_27, fixed_2_27, fixed_3_27, fixed_4_27, fixed_5_27])
# fixed_29 = np.array([fixed_1_29, fixed_2_29, fixed_3_29, fixed_4_29, fixed_5_29])
# fixed_31 = np.array([fixed_1_31, fixed_2_31, fixed_3_31, fixed_4_31, fixed_5_31])
fixed_42 = np.array([fixed_1_42, fixed_2_42, fixed_3_42, fixed_4_42, fixed_5_42])
# fixed_35 = np.array([fixed_1_35, fixed_2_35, fixed_3_35, fixed_4_35, fixed_5_35])
# fixed_57 = np.array([fixed_1_57, fixed_2_57, fixed_3_57, fixed_4_57, fixed_5_57])
# fixed_39 = np.array([fixed_1_39, fixed_2_39, fixed_3_39, fixed_4_39, fixed_5_39])
# fixed_45 = np.array([fixed_1_45, fixed_2_45, fixed_3_45, fixed_4_45, fixed_5_45])
# fixed_49 = np.array([fixed_1_49, fixed_2_49, fixed_3_49, fixed_4_49, fixed_5_49])
# fixed_53 = np.array([fixed_1_53, fixed_2_53, fixed_3_53, fixed_4_53, fixed_5_53])
# fixed_57 = np.array([fixed_1_57, fixed_2_57, fixed_3_57, fixed_4_57, fixed_5_57])
########################################################################################################################
# np.savetxt('stat_data_27.txt', fixed_27, delimiter=' ', header='Station data for year 2027')
# np.savetxt('stat_data_29.txt', fixed_29, delimiter=' ', header='Station data for year 2029')
# np.savetxt('stat_data_31.txt', fixed_31, delimiter=' ', header='Station data for year 2031')
np.savetxt('stat_data_42.txt', fixed_42, delimiter=' ', header='Station data for year 2042')
# np.savetxt('stat_data_35.txt', fixed_35, delimiter=' ', header='Station data for year 2035')
# np.savetxt('stat_data_57.txt', fixed_57, delimiter=' ', header='Station data for year 2057')
# np.savetxt('stat_data_39.txt', fixed_39, delimiter=' ', header='Station data for year 2039')
# np.savetxt('stat_data_45.txt', fixed_45, delimiter=' ', header='Station data for year 2045')
# np.savetxt('stat_data_49.txt', fixed_49, delimiter=' ', header='Station data for year 2049')
# np.savetxt('stat_data_53.txt', fixed_53, delimiter=' ', header='Station data for year 2053')
# np.savetxt('stat_data_57.txt', fixed_57, delimiter=' ', header='Station data for year 2057')
########################################################################################################################

