import matplotlib.pyplot as plt
import numpy as np
from numpy import *
# import cv2
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from dsw_dpt_slices import GRID, POR, dz, fluid_type, typo
from files_readme import celltopdpt, phi, cell_vol
from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
from files_readme import phi, cel_cen_dpt, dx, dy, dz, dx_accum, dy_accum, A, B

rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
k_range = 101

# Now we account for the different water saturations corresponding to the different times from year 2027 all the way
# to year 2039
################################################################################
################################################################################

# gv27_27 = np.loadtxt("grav_data_27.txt", delimiter=" ")


# The function calculating the gravity responses delta gz
def grav_mod(sat_1, sat_2):
    XSCALE = 84
    YSCALE = 89

    delta_g = np.zeros((B, A))
    fixed1, fixed2, fixed3, fixed4, fixed5 = None, None, None, None, None
    # rho_1 = ((1 - soil_1) * rho_water) + ((1 - sat_1) * rho_oil)
    # rho_2 = ((1 - soil_2) * rho_water) + ((1 - sat_2) * rho_oil)
    rho_1 = (sat_1 * rho_water) + ((1 - sat_1) * rho_oil)
    rho_2 = (sat_2 * rho_water) + ((1 - sat_2) * rho_oil)
    delta_rho = rho_2 - rho_1
    a = G * (cel_cen_dpt - 400) * delta_rho * phi
    b = dx * dy * dz  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    for y_step in range(0, B):
        for x_step in range(0, A):
            print("y", y_step)
            print("x", x_step)
            m = 0
            for k in range(k_range):
                for j in range(B):
                    for i in range(A):
                        r = np.sqrt(((x_step * XSCALE) - dx_accum[m]) ** 2 + ((y_step * YSCALE) - dy_accum[m]) ** 2)
                        d = pow((pow(r, 2) + pow(cel_cen_dpt[m] - 400, 2)), 1.5)
                        delta_g[158 - y_step, x_step] += ((c[m]) / d)
                        m += 1

            if x_step * XSCALE == 3948 and y_step * YSCALE == 11659:
                fixed1 = delta_g[157 - y_step, x_step]
            if x_step * XSCALE == 7812 and y_step * YSCALE == 6319:
                fixed2 = delta_g[157 - y_step, x_step]
            if x_step * XSCALE == 9324 and y_step * YSCALE == 8099:
                fixed3 = delta_g[157 - y_step, x_step]
            if x_step * XSCALE == 8652 and y_step * YSCALE == 4717:
                fixed4 = delta_g[157 - y_step, x_step]
            if x_step * XSCALE == 5796 and y_step * YSCALE == 7209:
                fixed5 = delta_g[157 - y_step, x_step]
    return delta_g * 1e8, fixed1 * 1e8, fixed2 * 1e8, fixed3 * 1e8, fixed4 * 1e8, fixed5 * 1e8
    # return delta_g * 1e8


################################################################################
print('Nooooooooooooooooooooooooooooooooooooooooooooooooooow')
# gv27_29 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_31 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_33 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_35 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_39 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_45 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_49 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_53 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_57 = grav_mod(sat27, sat29, soil_27, soil_29)

# gv27_27, fixed_1_27, fixed_2_27, fixed_3_27, fixed_4_27, fixed_5_27 = grav_mod(sat27, sat27, soil_27, soil_27)
# np.savetxt('grav_data_27_101_main.txt', gv27_27, delimiter=' ', header='Gravity data for year 2027')
# gv27_29, fixed_1_29, fixed_2_29, fixed_3_29, fixed_4_29, fixed_5_29 = grav_mod(sat27, sat29, soil_27, soil_29)
# np.savetxt('grav_data_29_101_main.txt', gv27_29, delimiter=' ', header='Gravity data for year 2029')
# gv27_31, fixed_1_31, fixed_2_31, fixed_3_31, fixed_4_31, fixed_5_31 = grav_mod(sat27, sat31, soil_27, soil_31)
# np.savetxt('grav_data_31_101_main.txt', gv27_31, delimiter=' ', header='Gravity data for year 2031')
# gv27_33, fixed_1_33, fixed_2_33, fixed_3_33, fixed_4_33, fixed_5_33 = grav_mod(sat27, sat33, soil_27, soil_33)
# np.savetxt('grav_data_33_101_main.txt', gv27_33, delimiter=' ', header='Gravity data for year 2033')
# gv27_35, fixed_1_35, fixed_2_35, fixed_3_35, fixed_4_35, fixed_5_35 = grav_mod(sat27, sat35, soil_27, soil_35)
# np.savetxt('grav_data_35_101_main.txt', gv27_35, delimiter=' ', header='Gravity data for year 2035')
# gv27_39, fixed_1_39, fixed_2_39, fixed_3_39, fixed_4_39, fixed_5_39 = grav_mod(sat27, sat39, soil_27, soil_39)
# np.savetxt('grav_data_39_101_main.txt', gv27_39, delimiter=' ', header='Gravity data for year 2039')
# gv27_45, fixed_1_45, fixed_2_45, fixed_3_45, fixed_4_45, fixed_5_45 = grav_mod(sat27, sat45, soil_27, soil_45)
# np.savetxt('grav_data_45_101_main.txt', gv27_45, delimiter=' ', header='Gravity data for year 2045')
# gv27_49, fixed_1_49, fixed_2_49, fixed_3_49, fixed_4_49, fixed_5_49 = grav_mod(sat27, sat49, soil_27, soil_49)
# np.savetxt('grav_data_49_101_main.txt', gv27_49, delimiter=' ', header='Gravity data for year 2049')
# gv27_53, fixed_1_53, fixed_2_53, fixed_3_53, fixed_4_53, fixed_5_53 = grav_mod(sat27, sat53, soil_27, soil_53)
# np.savetxt('grav_data_53_101_main.txt', gv27_53, delimiter=' ', header='Gravity data for year 2053')
gv27_57, fixed_1_57, fixed_2_57, fixed_3_57, fixed_4_57, fixed_5_57 = grav_mod(sat27, sat57)
np.savetxt('grav_data_57_101_main.txt', gv27_57, delimiter=' ', header='Gravity data for year 2057')
#
# fixed_27 = np.array([fixed_1_27, fixed_2_27, fixed_3_27, fixed_4_27, fixed_5_27])
# fixed_29 = np.array([fixed_1_29, fixed_2_29, fixed_3_29, fixed_4_29, fixed_5_29])
# fixed_31 = np.array([fixed_1_31, fixed_2_31, fixed_3_31, fixed_4_31, fixed_5_31])
# fixed_33 = np.array([fixed_1_33, fixed_2_33, fixed_3_33, fixed_4_33, fixed_5_33])
# fixed_35 = np.array([fixed_1_35, fixed_2_35, fixed_3_35, fixed_4_35, fixed_5_35])
# fixed_39 = np.array([fixed_1_39, fixed_2_39, fixed_3_39, fixed_4_39, fixed_5_39])
# fixed_45 = np.array([fixed_1_45, fixed_2_45, fixed_3_45, fixed_4_45, fixed_5_45])
# fixed_49 = np.array([fixed_1_49, fixed_2_49, fixed_3_49, fixed_4_49, fixed_5_49])
# fixed_53 = np.array([fixed_1_53, fixed_2_53, fixed_3_53, fixed_4_53, fixed_5_53])
# fixed_57 = np.array([fixed_1_57, fixed_2_57, fixed_3_57, fixed_4_57, fixed_5_57])
#
# np.savetxt('stat_data_27.txt', fixed_27, delimiter=' ', header='Station data for year 2027')
# np.savetxt('stat_data_29.txt', fixed_29, delimiter=' ', header='Station data for year 2029')
# np.savetxt('stat_data_31.txt', fixed_31, delimiter=' ', header='Station data for year 2031')
# np.savetxt('stat_data_33.txt', fixed_33, delimiter=' ', header='Station data for year 2033')
# np.savetxt('stat_data_35.txt', fixed_35, delimiter=' ', header='Station data for year 2035')
# np.savetxt('stat_data_39.txt', fixed_39, delimiter=' ', header='Station data for year 2039')
# np.savetxt('stat_data_45.txt', fixed_45, delimiter=' ', header='Station data for year 2045')
# np.savetxt('stat_data_49.txt', fixed_49, delimiter=' ', header='Station data for year 2049')
# np.savetxt('stat_data_53.txt', fixed_53, delimiter=' ', header='Station data for year 2053')
# np.savetxt('stat_data_57.txt', fixed_57, delimiter=' ', header='Station data for year 2057')
# #
# diff27_27 = gv27_27 - gv27_27
# diff27_29 = gv27_29 - gv27_27
# diff27_31 = gv27_31 - gv27_27
# diff27_33 = gv27_33 - gv27_27
# diff27_35 = gv27_35 - gv27_27
# diff27_39 = gv27_39 - gv27_27
# diff27_45 = gv27_45 - gv27_27
# diff27_49 = gv27_49 - gv27_27
# diff27_53 = gv27_53 - gv27_27
# diff27_57 = gv27_57 - gv27_27
# ################################################################################
# gv27_31 = grav_mod(sat27, sat31)
# gv27_33 = grav_mod(sat27, sat33)
# gv27_35 = grav_mod(sat27, sat35)
# gv27_45 = grav_mod(sat27, sat45)
# gv27_49 = grav_mod(sat27, sat49)
# gv27_53 = grav_mod(sat27, sat53)
# print("min27-57", np.amin(gv27_57))
# print("max27-57", np.amax(gv27_57))
# gv27 = np.reshape(gv27, (157, 159))
# print(gv27_57)################################################################################
# Plotting of the gravity modeling results
# fig = plt.figure()
# dx = linspace(-50, 50, num=100, endpoint=True)
# dy = linspace(-50, 50, num=100, endpoint=True)
# plt.plot(dx, dy)

# sns.set(style="white")
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111)
# nptsx, nptsy = 100, 100
# dxg, dyg = np.meshgrid(np.linspace(dx.min(), dx.max(), nptsx), np.linspace(dy.min(), dy.max(), nptsy))
# triangles = tri.Triangulation(dx, dy)
# tri_interp = tri.CubicTriInterpolator(triangles, gv27)
# gv27g = tri_interp(dxg, dyg)
# change levels here according to your data# fig = plt.figure()
# levels = np.linspace(0, 10, 5)
# colormap = ax.contourf(dxg, dyg, gv27, levels,
#                      cmap=plt.cm.Blues,
#                     norm=plt.Normalize(vmax=gv27.max(), vmin=gv27.min()))
# plot data points
# ax.plot(dx, dy, color="#444444", marker="o", linestyle="", markersize=10)

# add a colorbar
# fig.colorbar(colormap, orientation='vertical', shrink=0.85)  # horizontal colour bar

# graph extras: look at xlim and ylim
# ax.set_xlim((0, 10))
# ax.set_ylim((0, 10))
# ax.set_aspect("equal", "box")
# plt.show()

# plt.ylabel('$y$', fontsize=18)
# #fig.savefig('')
# plt.colorbar()
# color_map = plt.imshow(gv27)

# plt.style.use('classic')
# X = np.linspace(0, 13213.48, 157)
# Y = np.linspace(0, 14173.40, 159)
# XSCALE = 84
# YSCALE = 89
#
# dt_labels_tmstp = [gv27_27, gv27_29, gv27_31, gv27_33, gv27_35, gv27_39, gv27_45, gv27_49, gv27_53, gv27_57]
# # dt_labels_tmstp = [gv27_27, gv27_31, gv27_35, gv27_45, gv27_57]
# diff_dt_labels_tmstp = [diff27_27, diff27_29, diff27_31, diff27_33, diff27_35, diff27_39, diff27_45, diff27_49,
#                         diff27_53, diff27_57]
# labels_tmstp = ["gv_2727", "gv_2729", "gv_2731", "gv_2733", "gv_2735", "gv_2739", "gv_2745", "gv_2749",
#                 "gv_2753", "gv_2757"]

# diff_dt_labels_tmstp = [diff27_27, diff27_31, diff27_35, diff27_45, diff27_57]

# labels_tmstp = ["gv_2727", "gv_2731", "gv_2735", "gv_2745", "gv_2757"]

# for i in range(10):
#     dt = dt_labels_tmstp[i]
#     fig = plt.figure(figsize=(6.4, 4.8))
#     ax = fig.add_subplot(111)
#     ax.set_aspect('auto')
#     f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
#     plt.plot(3948 * 0.001, 11659 * 0.001, marker='s', markersize=7, color="w", linestyle='None', label="Station 1")
#     plt.plot(7812 * 0.001, 6319 * 0.001, marker='o', markersize=7, color="w", linestyle='None', label="Station 2")
#     plt.plot(9324 * 0.001, 8099 * 0.001, marker='^', markersize=7, color="w", linestyle='None', label="Station 3")
#     plt.plot(8652 * 0.001, 4717 * 0.001, marker='v', markersize=7, color="w", linestyle='None', label="Station 4")
#     plt.plot(5796 * 0.001, 7209 * 0.001, marker='>', markersize=7, color="w", linestyle='None', label="Station 5")
#     Z = f(Y * 0.001, X * 0.001)
#     Z = gaussian_filter(dt, sigma=0, mode='reflect')
#     plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
#                cmap=plt.cm.get_cmap('jet', 1000),
#                # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
#                interpolation='nearest', origin='upper')
#
#     clb = plt.colorbar()
#     tick_locator = ticker.MaxNLocator(nbins=7)
#     clb.locator = tick_locator
#     clb.update_ticks()
#     clb.ax.set_title('$(\u03BC$Gal)', fontsize=15)
#     ax.set_ylabel("Y(m)", labelpad=15)
#     plt.clim(0, 17)  #####################################################################################
#     # plt.title(labels_tmstp[i] + "_layerstack_" + phs)
#     plt.xlabel("X [$km$]", fontsize=18)
#     # ax.xaxis.tick_top()
#     plt.ylabel("Y [$km$]", fontsize=18)
#     ax = plt.gca()
#     legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
#     # ax.legend()
#     # Put a nicer background color on the legend.
#     # legend.get_frame()
#     # ax.set_xticks(X)
#     # ax.set_yticks(Y)
#     ax.grid(color='b', linestyle='-', linewidth=0.5)
#     # ax.set_xticklabels(np.arange(1, 158, 1))
#     # ax.set_yticklabels(np.arange(1, 160, 1))
#     plt.title(labels_tmstp[i])
#     plt.show()
#     fig.savefig(labels_tmstp[i])

# legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
# plt.legend(loc="upper right")
# plt.legend(loc=right, prop={'size': 6})
# Put a nicer background color on the legend.
# legend.get_frame()
########################################################################################################################
# diff = np.zeros(10)
# diff[0] = sum(gv27_27 - gv27_27)
# diff[1] = sum(gv27_29 - gv27_27)
# diff[2] = sum(gv27_31 - gv27_29)
# diff[3] = sum(gv27_33 - gv27_31)
# diff[4] = sum(gv27_35 - gv27_33)
# diff[5] = sum(gv27_39 - gv27_35)
# diff[6] = sum(gv27_45 - gv27_39)
# diff[7] = sum(gv27_49 - gv27_45)
# diff[8] = sum(gv27_53 - gv27_49)
# diff[9] = sum(gv27_57 - gv27_53)
# print(diff)

# time = [2027, 2029, 2031, 2033, 2035, 2039, 2045, 2049, 2053, 2057]
# time = [2027, 2031, 2035, 2045, 2057]
# time1 = [27, 29, 31, 33, 35, 39, 45, 49, 53, 57]

# Bar chart showing the gravity changes from one timestep to the next
# fig10, ax = plt.subplots()
# plt.grid()
# x = np.arange(1, len(diff) + 1)
# plt.bar(x, diff, align='center')
# plt.xticks(x, time1)
# plt.xlabel('$Timesteps (years)$', fontsize=13)
# plt.ylabel('$(\u03BCGal)$', fontsize=13)
# plt.show()
# fig10.savefig("overview.pdf")
#
# result1 = [fixed_1_27, fixed_1_29, fixed_1_31, fixed_1_33, fixed_1_35, fixed_1_39, fixed_1_45, fixed_1_49, fixed_1_53,
#            fixed_1_57]
# result2 = [fixed_2_27, fixed_2_29, fixed_2_31, fixed_2_33, fixed_2_35, fixed_2_39, fixed_2_45, fixed_2_49, fixed_2_53,
#            fixed_2_57]
# result3 = [fixed_3_27, fixed_3_29, fixed_3_31, fixed_3_33, fixed_3_35, fixed_3_39, fixed_3_45, fixed_3_49, fixed_3_53,
#            fixed_3_57]
# result4 = [fixed_4_27, fixed_4_29, fixed_4_31, fixed_4_33, fixed_4_35, fixed_4_39, fixed_4_45, fixed_4_49, fixed_4_53,
#            fixed_4_57]
# result5 = [fixed_5_27, fixed_5_29, fixed_5_31, fixed_5_33, fixed_5_35, fixed_5_39, fixed_5_45, fixed_5_49, fixed_5_53,
#            fixed_5_57]

# result1 = [fixed_1_27, fixed_1_31, fixed_1_35, fixed_1_45, fixed_1_57]
# result2 = [fixed_2_27, fixed_2_31, fixed_2_35, fixed_2_45, fixed_2_57]
# result3 = [fixed_3_27, fixed_3_31, fixed_3_35, fixed_3_45, fixed_3_57]

# fig11, ax = plt.subplots()
# x = np.arange(10)
# ax.plot(result1, 'k', label='Station 1')
# ax.plot(result2, 'k--', label='Station 2')
# ax.plot(result3, 'k:', label='Station 3')
# ax.plot(result4, 'k-.', label='Station 4')
# ax.plot(result5, '^k', label='Station 5')
# legend = ax.legend(loc='upper right', shadow=True, fontsize='small')
# # Put a nicer background color on the legend.
# legend.get_frame()
# # plt.title('Gravity changes for three measuring stations')
# plt.xlabel('$Timesteps (years)$', fontsize=18)
# plt.ylabel('$(\u03BCGal)$', fontsize=18)
# plt.xticks(x, time1)
# plt.show()
# fig11.savefig("final.pdf")
#
# print('Result1', result1)
# print('Result2', result2)
# print('Result3', result3)
# print('Result4', result4)
# print('Result5', result5)
