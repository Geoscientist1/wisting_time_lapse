import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from gravity_inversion_Wisting import M, x_range, y_range, z_range, x_range_data, y_range_data
from Grav_Modeling_coarse_grid import *

# from dsw_dpt_slices import GRID, POR, dz, fluid_type, typo
# from files_readme import celltopdpt, phi, cell_vol
# from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
# from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
# from files_readme import phi, cel_cen_dpt, dx, dy, dz, dx_accum, dy_accum, A, B

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
# gv27_29, fixed_1_29, fixed_2_29, fixed_3_29, fixed_4_29, fixed_5_29 = grav_mod(sat27, sat29, soil_27, soil_29)
# gv27_31, fixed_1_31, fixed_2_31, fixed_3_31, fixed_4_31, fixed_5_31 = grav_mod(sat27, sat31, soil_27, soil_31)
# gv27_33, fixed_1_33, fixed_2_33, fixed_3_33, fixed_4_33, fixed_5_33 = grav_mod(sat27, sat33, soil_27, soil_33)
# gv27_35, fixed_1_35, fixed_2_35, fixed_3_35, fixed_4_35, fixed_5_35 = grav_mod(sat27, sat35, soil_27, soil_35)
# gv27_39, fixed_1_39, fixed_2_39, fixed_3_39, fixed_4_39, fixed_5_39 = grav_mod(sat27, sat39, soil_27, soil_39)
# gv27_45, fixed_1_45, fixed_2_45, fixed_3_45, fixed_4_45, fixed_5_45 = grav_mod(sat27, sat45, soil_27, soil_45)
# gv27_49, fixed_1_49, fixed_2_49, fixed_3_49, fixed_4_49, fixed_5_49 = grav_mod(sat27, sat49, soil_27, soil_49)
# gv27_53, fixed_1_53, fixed_2_53, fixed_3_53, fixed_4_53, fixed_5_53 = grav_mod(sat27, sat53, soil_27, soil_53)
# gv27_57 = grav_mod(sat27, sat57, soil_27, soil_57)


# gv27_57 = np.loadtxt("grav_data_57_101_coarse_dt_point_grid_5.txt", delimiter=" ")
#
# plt.style.use('classic')
# X = np.linspace(0, 13213.48, 157)
# Y = np.linspace(0, 14173.40, 159)
# XSCALE = 84
# YSCALE = 89
#
# dt_labels_tmstp = [gv27_57]
# labels_tmstp = ["gv_2757"]
#
# dt = dt_labels_tmstp[0]
# fig = plt.figure(figsize=(6.4, 4.8))
# ax = fig.add_subplot(111)
# ax.set_aspect('auto')
# f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
# Z = f(Y * 0.001, X * 0.001)
# Z = gaussian_filter(dt, sigma=0, mode='reflect')
# plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
#            cmap=plt.cm.get_cmap('jet', 1000),
#            # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
#            interpolation='nearest', origin='upper')
#
# clb = plt.colorbar()
# tick_locator = ticker.MaxNLocator(nbins=7)
# clb.locator = tick_locator
# clb.update_ticks()
# clb.ax.set_title('$(\u03BC$Gal)', fontsize=15)
# ax.set_ylabel("Y(m)", labelpad=15)
# plt.clim(0, 100)
# # plt.title(labels_tmstp[i] + "_layerstack_" + phs)
# plt.title("Non-Uniform Grid")
# plt.xlabel("X [$km$]")
# # ax.xaxis.tick_top()
# plt.ylabel("Y [$km$]")
# ax = plt.gca()
# ax.patch.set_visible(False)
# fig.patch.set_visible(False)
# legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
# ax.legend()
# Put a nicer background color on the legend.
# legend.get_frame()
# ax.set_xticks(X)
# ax.set_yticks(Y)
# ax.grid(color='b', linestyle='-', linewidth=0.5)
# ax.set_xticklabels(np.arange(1, 158, 1))
# ax.set_yticklabels(np.arange(1, 160, 1))
# plt.show()
# fig.savefig(labels_tmstp[i], transparent=True)

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
# np.savetxt('grav_data_27_101.txt', gv27_27, delimiter=' ', header='Gravity data for year 2027')
# np.savetxt('grav_data_29_101.txt', gv27_29, delimiter=' ', header='Gravity data for year 2029')
# np.savetxt('grav_data_31_101.txt', gv27_31, delimiter=' ', header='Gravity data for year 2031')
# np.savetxt('grav_data_33_101.txt', gv27_33, delimiter=' ', header='Gravity data for year 2033')
# np.savetxt('grav_data_35_101.txt', gv27_35, delimiter=' ', header='Gravity data for year 2035')
# np.savetxt('grav_data_39_101.txt', gv27_39, delimiter=' ', header='Gravity data for year 2039')
# np.savetxt('grav_data_45_101.txt', gv27_45, delimiter=' ', header='Gravity data for year 2045')
# np.savetxt('grav_data_49_101.txt', gv27_49, delimiter=' ', header='Gravity data for year 2049')
# np.savetxt('grav_data_53_101.txt', gv27_53, delimiter=' ', header='Gravity data for year 2053')
# np.savetxt('grav_data_57_101.txt', gv27_57, delimiter=' ', header='Gravity data for year 2057')
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

# import matplotlib.pyplot as plt
# import numpy as np
# from numpy import *
#
# A = 3
# B = 5
# k_range = 7
# k_range_new = 2
# temp = 0
#
# arr = np.ones((k_range, B, A))
# arr1 = np.zeros((k_range_new, B, A))
# print("arr is: ", arr)
# print("arr shape is:", arr.shape)
# for i in range(A):
#     for j in range(B):
#         s = 0
#         for k in range(k_range):
#             temp += arr[k, j, i]
#             if k != 0 and mod(k, 3) == 0:
#                 arr1[s, j, i] = temp / 7
#                 temp = 0
#                 s += 1
#
# print("This is arr: ", arr)
# print("This is arr1: ", arr1)
#
# arr1 = arr1.ravel()
# print("arr1 is: ", arr1)
# arr2 = np.reshape(arr1, (k_range_new, B, A))
# print("arr2 is:", arr2)

# if arr == arr2.all():
#     print("The two arrays are the same")
# else:
#     print("The two arrays are not the same")
# print("arr1 shape is:", arr1.shape)


#
# # from scipy.ndimage import gaussian_filter
# # from scipy import misc
# # import matplotlib.pyplot as plt
# # from files_readme import sat27, sat29, poro, porv, cell_vol, dx, dy, dz, vg_57, vo_27, vo_26, vo_45, vw_45, vw_26, \
# #     vw_57, \
# #     sat27, \
# #     vo_57
# # #from Grav_rev import val27_57, val27_29
#
# # rho_oil = 840  # Oil density
# # rho_water = 1040  # Water density
# # rho_1 = ((1 - sat27) * rho_oil) + (sat27 * rho_water)
# # rho_2 = ((1 - sat29) * rho_oil) + (sat29 * rho_water)
# # delta_rho = rho_2 - rho_1
# # S = np.zeros(len(delta_rho))
# # R = np.zeros((157, 159))
# # for i in range(157):
# #     for j in range(159):
# #         for k in range(101):
# #             m = i * j * k
# #             print(delta_rho[m])
# #             S[m] = delta_rho[m]
# #             R[i, j] += S[m]
# #
# # # delta_rho = np.reshape(delta_rho, (157, 159))
# # np.savetxt('RHO.txt', S)
# # np.savetxt('BHO.txt', R)
#
# print("##########################################################")
#
# # plt.style.use('classic')
# # fig = plt.figure(figsize=(6.4, 4.8))
# # ax = fig.add_subplot(111)
# # ax.set_aspect('equal')
# # plt.imshow(val27_57, cmap=plt.cm.get_cmap('jet', 75), interpolation='nearest', origin='lower')
# # plt.colorbar(extend='max')
# # plt.clim(0, 3000)
# # plt.title('rho_27_57')
# # plt.xlabel("X(m)")
# # plt.ylabel("Y(m)")
# # plt.show()
# # fig.savefig("rho27_27.png")
# #
# # fig2 = plt.figure(figsize=(6.4, 4.8))
# # ax = fig2.add_subplot(111)
# # ax.set_aspect('equal')
# # result2 = gaussian_filter(val27_57, sigma=1)
# # plt.imshow(result2, origin='lower')
# # plt.colorbar(extend='max')
# # plt.clim(0, 2800)
# # plt.show()
# #
# print("##########################################################")
#
# # disx = 0
# # disy = 0
# # disz = 0
# # i = 0
# # j = 0
# #
# # for n in range(2521263):
# #     if dy[n] != 0 and i < 159:
# #         disy += dy[n]
# #         i += 1
# #         if disy > 14316 or i > 159:
# #             break
# #
# # for n in range(2521263):
# #     if dz[n] != 0:
# #         disz += dz[n]
# #         j += 1
# #         if disz > 487:
# #             break
# #
# # print("Total distance dy: ", disy, i)
# # print("Total distance dy: ", disz, j)
#
# ########################################################################################################################
# # def avg(list1, nx, ny, nz):
# #     temp1 = 0
# #     temp2 = 0
# #     N = len(list1)
# #     n = 0
# #     new = []
# #     new1 = []
# #     m = 0
# #     for i in range(nx):
# #         for j in range(ny):
# #             for k in range(nz):
# #                 if dx[m] == dy[m] == dz[m] == 0:
# #                     temp1 = 0
# #                 else:
# #                     temp1 += list1[m]
# #                     if np.mod(k, 13) == 0:
# #                         new.append(temp1 / 14)
# #                         temp1 = 0
# #                 m += 1
# #     print(len(new))
# #     s = 0
# #     for i in range(nx):
# #         for j in range(ny):
# #             for k in range(12):
# #                 temp2 += new[s]
# #                 if k == 11:
# #                     new1.append(temp2 / 12)
# #                     temp2 = 0
# #                 s += 1
# #     return np.array(new1)
# #
# #
# # print(avg(dx, 157, 159, 101))
#
# ########################################################################################################################
#
#
# # # Read GRDECL File
# # Model = GeologyModel(filename='./ExampleData/dome.grdecl')
# #
# # # Convert ECLIPSE grdecl format into VTK
# # Model.GRDECL2VTK()
# #
# # # Decompose the model into sub-volumes in terms of faults automatically (this function requires shapely library)
# # Model.decomposeModel()
# #
# # # Output to VTK format
# # Model.Write2VTU()
# #
# # # Load a custom new keyword from file
# # TempData = Model.LoadCellData(varname="TEMP", filename='./ExampleData/dome_Temperature.txt')
# #
# # # Update model and output to VTK format
# # Model.Update()
# # Model.Write2VTU()
#
# #
# # a = 0.81  # Turtuosity factor in moderately porous sand
# # rw = 0.05  # Formation water resistivity
# # x = 2.1  # Saturation exponent in the oil zone
# # y = 2  # Cementation factor
# #
# #
# # def Sw_profile(sat_1, dz):
# #     Sw = list()
# #     tvd = list()
# #     cell_x = 114
# #     cell_y = 59
# #     m = 0
# #     for k in range(101):
# #         for j in range(159):
# #             for i in range(157):
# #                 if i == cell_x and j == cell_y and dz[m] != 0:
# #                     Sw.append(sat_1[m])
# #                     tvd.append(dz[m])
# #                 m += 1
# #
# #     Sw = np.array(Sw, dtype=object)
# #     tvd = np.array(tvd, dtype=object)
# #     tvd = np.cumsum(tvd) + 693.65
# #     return Sw, tvd
# #
# #
# # def Rt_profile(sat_1, dz):
# #     tvd = list()
# #     rt = list()
# #     porosity = list()
# #     cell_x = 114
# #     cell_y = 59
# #     m = 0
# #     for k in range(101):
# #         for j in range(159):
# #             for i in range(157):
# #                 if i == cell_x and j == cell_y and dz[m] != 0:
# #                     tvd.append(dz[m])
# #                     rt.append((a * rw) * (1 / (power(sat_1[m], x) * power(phi[m], y))))
# #                     porosity.append(phi[m])
# #                 m += 1
# #
# #     rt = np.array(rt, dtype=object)
# #     porosity = np.array(porosity, dtype=object)
# #     tvd = np.array(tvd, dtype=object)
# #     tvd = np.cumsum(tvd) + 693.65
# #     return rt, tvd, porosity
# #
# #
# # Sw_27, TVD = Sw_profile(sat27, dz)
# # Rt_27, TVD, porosity = Rt_profile(sat27, dz)
# #
# # fig = plt.figure(figsize=(19.05/2.54, 23.01/2.54), facecolor='#e1ddbf')
# # ax = fig.add_subplot(111)
# # plt.plot(Sw_27, TVD, label='Sw', color='tab:cyan', linewidth=3, marker='h', markerfacecolor='lightgreen',
# #          markeredgewidth=1,
# #          markersize=7, markevery=3)
# # plt.plot(porosity, TVD, label='Porosity', color='tab:blue', linewidth=3, marker='h', markerfacecolor='lightskyblue',
# #          markeredgewidth=1,
# #          markersize=7, markevery=3)
# # plt.ylim((650, 800))
# # plt.gca().invert_yaxis()
# # plt.title('Water saturation along Well 7324/8-1')
# # plt.xlabel("Sw (Water Saturation)", fontsize=20)
# # plt.ylabel("TVD [m]", fontsize=20)
# # plt.setp(ax.spines.values(), linewidth=3)
# # # The ticks
# # ax.xaxis.set_tick_params(width=3)
# # ax.yaxis.set_tick_params(width=3)
# # plt.locator_params(axis='x')
# # plt.axvline(x=np.amin(Sw_27), color='r', linewidth=3, label='min Sw =' + str(np.amin(Sw_27)))
# # plt.axvline(x=np.amin(porosity), color='y', linewidth=3, label='min porosity =' + str(np.amin(porosity)))
# # plt.axhline(y=662, color='b', linewidth=3, label='Top Sto')
# # plt.legend(loc='best')
# # plt.show()
# # fig.savefig('Water_saturation_Profile')
# #
# #
# # fig = plt.figure(figsize=(19.05/2.54, 23.01/2.54), facecolor='#e1ddbf')
# # ax = fig.add_subplot(111)
# # plt.plot(Rt_27, TVD, label='Rt', color='k', linewidth=3, marker='h', markerfacecolor='lightgreen',
# #          markeredgewidth=1,
# #          markersize=7, markevery=3)
# # plt.ylim((650, 800))
# # plt.gca().invert_yaxis()
# # plt.title('Deep Resistivity along Well 7324/8-1')
# # plt.xlabel("Rt (Deep Resistivity)", fontsize=20)
# # plt.ylabel("TVD [m]", fontsize=20)
# # plt.xscale(value='log')
# # plt.setp(ax.spines.values(), linewidth=3)
# # # The ticks
# # ax.xaxis.set_tick_params(width=3)
# # ax.yaxis.set_tick_params(width=3)
# # ax.locator_params(axis='x')
# # plt.axvline(x=np.amax(Rt_27), color='r', linewidth=3, label='max Rt =' + str(np.amax(Rt_27)))
# # plt.axhline(y=662, color='b', linewidth=3, label='Top Sto')
# # plt.legend(loc='best')
# # plt.show()
# # fig.savefig('Deep_Resistivity_Profile')
#
#
# # def petrel2sblwiz(file, new_file):
# #     line_count = 0
# #     imported_model = np.zeros((101, 159, 157))
# #     bck_res = np.zeros((101, 159, 157))
# #
# #     with open(file, "r") as file:
# #         lines = file.readlines()
# #         file.close()
# #
# #     del lines[0]
# #     del lines[0]
# #     del lines[0]
# #     del lines[0]
# #     del lines[0]
# #     del lines[0]
# #     del lines[0]
# #     del lines[0]
# #     del lines[0]
# #
# #     new_file = open(new_file, "w+")
# #
# #     for line in lines:
# #         new_file.write(line)
# #
# #     new_file.close()
# #
# #     return None
# #
# #
# # resis_target = 'moon.txt'
# # resis_target_edited = 'moon_edited.txt'
# # petrel2sblwiz(resis_target, resis_target_edited)
#
#
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.tri as mtri
#
# # The values ​​related to each point. This can be a "Dataframe pandas"
# # for example where each column is linked to a variable <-> 1 dimension.
# # The idea is that each line = 1 pt in 4D.
# do_random_pt_example = True
#
# index_x = 0
# index_y = 1
# index_z = 2
# index_c = 3
# list_name_variables = ['x', 'y', 'z', 'c']
# name_color_map = 'jet'
#
# if do_random_pt_example:
#     number_of_points = 200
#     x = np.random.rand(number_of_points)
#     y = np.random.rand(number_of_points)
#     z = np.random.rand(number_of_points)
#     c = np.random.rand(number_of_points)
# else:
#     # Example where we have a "Pandas Dataframe" where each line = 1 pt in 4D.
#     # We assume here that the "data frame" "df" has already been loaded before.
#     x = df[list_name_variables[index_x]]
#     y = df[list_name_variables[index_y]]
#     z = df[list_name_variables[index_z]]
#     c = df[list_name_variables[index_c]]
# # end
# # -----
#
# # We create triangles that join 3 pt at a time and where their colors will be
# # determined by the values ​​of their 4th dimension. Each triangle contains 3
# # indexes corresponding to the line number of the points to be grouped.
# # Therefore, different methods can be used to define the value that
# # will represent the 3 grouped points and I put some examples.
# triangles = mtri.Triangulation(x, y).triangles
#
# choice_calcuation_colors = 1
# if choice_calcuation_colors == 1:  # Mean of the "c" values of the 3 pt of the triangle
#     colors = np.mean([c[triangles[:, 0]], c[triangles[:, 1]], c[triangles[:, 2]]], axis=0)
# elif choice_calcuation_colors == 2:  # Mediane of the "c" values of the 3 pt of the triangle
#     colors = np.median([c[triangles[:, 0]], c[triangles[:, 1]], c[triangles[:, 2]]], axis=0)
# elif choice_calcuation_colors == 3:  # Max of the "c" values of the 3 pt of the triangle
#     colors = np.max([c[triangles[:, 0]], c[triangles[:, 1]], c[triangles[:, 2]]], axis=0)
# # end
# # ----------
# # Displays the 4D graphic.
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# triang = mtri.Triangulation(x, y, triangles)
# surf = ax.plot_trisurf(triang, z, cmap=name_color_map, shade=False, linewidth=0.2)
# surf.set_array(colors)
# surf.autoscale()
#
# # Add a color bar with a title to explain which variable is represented by the color.
# cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
# cbar.ax.get_yaxis().labelpad = 15
# cbar.ax.set_ylabel(list_name_variables[index_c], rotation=270)
#
# # Add titles to the axes and a title in the figure.
# ax.set_xlabel(list_name_variables[index_x])
# ax.set_ylabel(list_name_variables[index_y])
# ax.set_zlabel(list_name_variables[index_z])
# plt.title('%s in function of %s, %s and %s' % (
#     list_name_variables[index_c], list_name_variables[index_x], list_name_variables[index_y],
#     list_name_variables[index_z]))
#
# plt.show()


# Import data
# import time
# import numpy as np
#
# from skimage import io
#


dt_diff = np.subtract(sat57, sat27)
p = []
for i in range(M):
    if dt_diff[i] != 0:
        p.append(dt_diff[i])
p = np.array(p)

map_mtr = np.zeros((M, len(p)))
j = 0
for i in range(M):
    if dt_diff[i] != 0:
        map_mtr[i, j] = 1
        j += 1

print("Alhamdulilah")

new_m = map_mtr @ p
gv1 = grav_mod_coarse(new_m)
np.savetxt('grav_data_coarse_grid_mapping_matrix_test.txt', gv1, delimiter=' ', header='Gravity data')
gv0 = np.loadtxt('grav_data_coarse_grid_coarse_dtpoint_grid.txt', delimiter=' ')
gv_diff = np.subtract(gv1, gv0)
np.savetxt('grav_data_difference.txt', gv_diff, delimiter=' ', header='Gravity data')

result = np.isin(dt_diff, new_m)
print(result)
jcb = np.loadtxt("jacobian_perturbation_method.txt", delimiter=" ")
#
#
#
# plt.style.use('classic')
#
# dt_labels_tmstp = [jcb]
# labels_tmstp = ["jacobian"]
#
# dt = dt_labels_tmstp[0]
# fig = plt.figure(figsize=(6.4, 4.8))
# ax = fig.add_subplot(111)
# ax.set_aspect('auto')
# plt.imshow(dt, cmap=plt.cm.get_cmap('jet', 1000),
#            # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
#            interpolation='nearest', origin='upper')
#
# clb = plt.colorbar()
# tick_locator = ticker.MaxNLocator(nbins=7)
# clb.locator = tick_locator
# clb.update_ticks()
# clb.ax.set_title('', fontsize=15)
# ax.set_ylabel("Y(m)", labelpad=15)
# plt.clim(0, 20)
# plt.axis('off')
# plt.title(labels_tmstp[i] + "_layerstack_" + phs)
# plt.title("")
# plt.xlabel("")
# ax.xaxis.tick_top()
# plt.ylabel("")
# ax = plt.gca()
# ax.patch.set_visible(False)
# fig.patch.set_visible(False)
# legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
# ax.legend()
# Put a nicer background color on the legend.
# legend.get_frame()
# ax.set_xticks(X)
# ax.set_yticks(Y)
# ax.grid(color='b', linestyle='-', linewidth=0.5)
# ax.set_xticklabels(np.arange(1, 158, 1))
# ax.set_yticklabels(np.arange(1, 160, 1))
# plt.show()
# fig.savefig(labels_tmstp[i], transparent=True)



# mesh = CreateTensorMesh(
#     origin=[793000, 9192500, 2690],
#     xcellstr="1000 500 50*250 500 1000",
#     ycellstr="1000 500 55*250 500 1000",
#     zcellstr="30*100.0 5*250.0 500",
# ).apply()
#
# mesh.plot(show_grid=True, color=True, show_edges=True)
########################################################################################################################

# reader.SetFilename("sat57_start_model_Wisting.txt")
# reader.update()
# my_vtk_dataset = reader.GetOutput()
# mesh1 = CreateTensorMesh(
#     origin=[0, 0, 0],
#     xcellstr="157*87",
#     ycellstr="159*89",
#     zcellstr="5*70",
#     data_name=dt,
# ).apply()

# mesh1.plot(show_grid=True, color=True, show_edges=True)
# output = XYZTextReader().apply('sat57_start_model_Wisting.txt')
# centers = AppendCellCenters().apply(output)
# # centers
# centers.plot()

import numpy as np
import time
import plotly.graph_objects as go
from plotly.offline import plot  # for IDE use

# from paraview import util
#
# util.SetOutputWholeExtent(self, [0, 157, 0, 159, 0, 5])
# a = np.loadtxt("sat57_start_model_Wisting.txt").flatten()
#
# output.SetDimensions(157, 159, 5)
# output.SetOrigin(0.0, 0.0, 0.0)
# output.SetSpacing(1, 1, 1)
# output.PointData.append(a, 'height')
x_range = 157  # This is the # of grid cells in the x-direction
y_range = 159  # This is the # of grid cells in the y-direction
z_range = 5
x_max = 14000
y_max = 14000
z_max = -1000


def multiple3dsurfaceplots(model):
    model_reshaped = np.reshape(model, (z_range, y_range, x_range))
    X = np.linspace(0, x_max, x_range)
    Y = np.linspace(0, y_max, y_range)
    # Z = np.linspace(z_max, -300, z_range)
    Z = np.linspace(-600, z_max, z_range)
    model_to_be_drawn = np.zeros((z_range, y_range, x_range))
    x = np.zeros((z_range, y_range, x_range))
    y = np.zeros((z_range, y_range, x_range))
    z = np.zeros((z_range, y_range, x_range))

    # 20 cells for the water layer with uniform resistivity being equal to 0.1 ohm-m
    for k in range(0, z_range):
        for j in range(y_range):
            for i in range(x_range):
                model_to_be_drawn[k, j, i] = model_reshaped[k, j, i]
                z[k, j, i] = Z[k]
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]

    return model_to_be_drawn, x, y, z


sat_model = np.loadtxt("sat57_coarse.txt", delimiter=" ")

fake_model, x_coordinates, y_coordinates, z_coordinates = multiple3dsurfaceplots(sat_model)

# ame = 'eye = (x:2, y:2, z:0.1)'
# camera = dict(
#         eye=dict(x=2, y=2, z=0.1)
#     )
# density_models = [fake_model]
# plot_title_labels = ["Test Model"]
# for i in range(2):
#     time.sleep(1)
#     dt = density_models[i]
#     dt_Smooth = gaussian_filter(dt, sigma=0, mode='reflect')
#     fig = go.Figure(data=[
#     go.Volume(x=x_coordinates.flatten(), y=y_coordinates.flatten(), z=z_coordinates.flatten(),
#                 value=dt_Smooth.flatten(), cmin=0, cmax=1, colorscale='Jet', opacity=.5, colorbar_len=0.3),
#         ])
#     fig.update_layout(scene_camera=camera, title=plot_title_labels[i], width=700, height=700)
#     plot(fig, auto_open=True)
#     fig.write_image(labels_tmstp_0[i] + ".png")
#     layout = go.Layout(
#         )



