import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from matplotlib import rcParams

########################################################################################################################

# This script plots the modeled ATR maps which are outputted by atr_gen.py
# The modeled atr maps are saved and hence running and computional time is saved by this script

########################################################################################################################
# def grav_data_plot(gravdata):
#     M = 157
#     N = 159
#     line_count = 0
#     output = np.zeros((N, M))
#
#     with open(gravdata, "r") as file:
#         for line in file:
#             if line != "\n":
#                 line_count += 1
#
#             arr = line.split()
#             for k in range(len(arr)):
#                 i = int(arr[3])
#                 j = 318 - int(arr[4])
#                 resist[j, i] = float(arr[2])
#
#     return resist, bck_res


gv27_27 = np.loadtxt("grav_data_27_101_main.txt", delimiter=" ")
gv27_29 = np.loadtxt("grav_data_29_101_main.txt", delimiter=" ")
gv27_31 = np.loadtxt("grav_data_31_101_main.txt", delimiter=" ")
gv27_32 = np.loadtxt("grav_data_32_101_main.txt", delimiter=" ")
gv27_33 = np.loadtxt("grav_data_33_101_main.txt", delimiter=" ")
gv27_35 = np.loadtxt("grav_data_35_101_main.txt", delimiter=" ")
gv27_37 = np.loadtxt("grav_data_37_101_main.txt", delimiter=" ")
gv27_39 = np.loadtxt("grav_data_39_101_main.txt", delimiter=" ")
gv27_41 = np.loadtxt("grav_data_41_101_main.txt", delimiter=" ")
gv27_42 = np.loadtxt("grav_data_42_101_main.txt", delimiter=" ")
gv27_43 = np.loadtxt("grav_data_43_101_main.txt", delimiter=" ")
gv27_45 = np.loadtxt("grav_data_45_101_main.txt", delimiter=" ")
gv27_47 = np.loadtxt("grav_data_47_101_main.txt", delimiter=" ")
gv27_49 = np.loadtxt("grav_data_49_101_main.txt", delimiter=" ")
gv27_51 = np.loadtxt("grav_data_51_101_main.txt", delimiter=" ")
gv27_53 = np.loadtxt("grav_data_53_101_main.txt", delimiter=" ")
gv27_55 = np.loadtxt("grav_data_55_101_main.txt", delimiter=" ")
gv27_57 = np.loadtxt("grav_data_57_101_main.txt", delimiter=" ")

result1 = [gv27_27[34, 42], gv27_29[34, 42], gv27_31[34, 42], gv27_33[34, 42], gv27_35[34, 42], gv27_37[34, 42],
           gv27_39[34, 42], gv27_41[34, 42], gv27_43[34, 42], gv27_45[34, 42], gv27_47[34, 42], gv27_49[34, 42],
           gv27_51[34, 42], gv27_53[34, 42], gv27_55[34, 42], gv27_57[34, 42]]

result2 = [gv27_27[54, 100], gv27_29[54, 100], gv27_31[54, 100], gv27_33[54, 100], gv27_35[54, 100],
           gv27_37[54, 100], gv27_39[54, 100], gv27_41[54, 100], gv27_43[54, 100], gv27_45[54, 100],
           gv27_47[54, 100], gv27_49[54, 100], gv27_51[54, 100], gv27_53[54, 100], gv27_55[54, 100],
           gv27_57[54, 100]]

result3 = [gv27_27[78, 109], gv27_29[78, 109], gv27_31[78, 109], gv27_33[78, 109], gv27_35[78, 109],
           gv27_37[78, 109], gv27_39[78, 109], gv27_41[78, 109], gv27_43[78, 109], gv27_45[78, 109],
           gv27_47[78, 109], gv27_49[78, 109], gv27_51[78, 109], gv27_53[78, 109], gv27_55[78, 109],
           gv27_57[78, 109]]

result4 = [gv27_27[51, 134], gv27_29[51, 134], gv27_31[51, 134], gv27_33[51, 134], gv27_35[51, 134],
           gv27_37[51, 134], gv27_39[51, 134], gv27_41[51, 134], gv27_43[51, 134], gv27_45[51, 134],
           gv27_47[51, 134], gv27_49[51, 134], gv27_51[51, 134], gv27_53[51, 134], gv27_55[51, 134],
           gv27_57[51, 134]]

result5 = [gv27_27[97, 147], gv27_29[97, 147], gv27_31[97, 147], gv27_33[97, 147], gv27_35[97, 147],
           gv27_37[97, 147], gv27_39[97, 147], gv27_41[97, 147], gv27_43[97, 147], gv27_45[97, 147],
           gv27_47[97, 147], gv27_49[97, 147], gv27_51[97, 147], gv27_53[97, 147], gv27_55[97, 147],
           gv27_57[97, 147]]

plt.style.use('classic')
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
XSCALE = 84
YSCALE = 89

dt_labels_tmstp = [gv27_32, gv27_42, gv27_57]

labels_tmstp = ["2032", "2042", "2057"]

# dt_labels_tmstp = [gv27_27, gv27_29, gv27_31, gv27_33, gv27_35, gv27_37, gv27_39, gv27_41, gv27_43, gv27_45, gv27_47,
#                    gv27_49, gv27_51, gv27_53, gv27_55, gv27_57]
# labels_tmstp = ["gv_2727", "gv_2729", "gv_2731", "gv_2733", "gv_2735", "gv_2737", "gv_2739", "gv_2741", "gv_2743",
#                 "gv_2745", "gv_2747", "gv_2749", "gv_2751", "gv_2753", "gv_2755", "gv_2757"]


for i in range(3):
    rcParams['font.weight'] = 'bold'
    dt = dt_labels_tmstp[i]
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')
    f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
    plt.plot(3512 * 0.001, 11103 * 0.001, marker='s', markersize=7, color="w", linestyle='None', label="Station 1")
    plt.plot(8400 * 0.001, 9345 * 0.001, marker='>', markersize=7, color="w", linestyle='None', label="Station 2")
    plt.plot(9156 * 0.001, 7209 * 0.001, marker='v', markersize=7, color="w", linestyle='None', label="Station 3")
    plt.plot(11256 * 0.001, 9612 * 0.001, marker='^', markersize=7, color="w", linestyle='None', label="Station 4")
    plt.plot(12348 * 0.001, 5518 * 0.001, marker='o', markersize=7, color="w", linestyle='None', label="Station 5")
    Z = f(Y * 0.001, X * 0.001)
    Z = gaussian_filter(dt, sigma=0, mode='reflect')
    plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
               cmap=plt.cm.get_cmap('jet', 1000),
               # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
               interpolation='nearest', origin='upper')

    clb = plt.colorbar()
    tick_locator = ticker.MaxNLocator(nbins=7)
    clb.locator = tick_locator
    clb.update_ticks()
    clb.ax.set_title('$(\u03BC$Gal)', fontsize=15, weight='bold')
    ax.set_ylabel("Y(m)", labelpad=15)
    plt.clim(0, 13)
    # plt.title(labels_tmstp[i] + "_layerstack_" + phs)
    plt.xlabel("East [km]", fontsize=15, weight='bold')
    # ax.xaxis.tick_top()
    plt.ylabel("North [km]", fontsize=15, weight='bold')
    ax = plt.gca()
    ax.patch.set_visible(False)
    fig.patch.set_visible(False)
    legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
    # ax.legend()
    # Put a nicer background color on the legend.
    # legend.get_frame()
    # ax.set_xticks(X)
    # ax.set_yticks(Y)
    # ax.grid(color='b', linestyle='-', linewidth=0.5)
    # ax.set_xticklabels(np.arange(1, 158, 1))
    # ax.set_yticklabels(np.arange(1, 160, 1))
    plt.title(labels_tmstp[i], fontsize=15, weight='bold')
    plt.show()
    # fig.savefig(labels_tmstp[i], transparent=True)

    diff = np.zeros(16)
    diff[0] = sum(gv27_27 - gv27_27)
    diff[1] = sum(gv27_29 - gv27_27)
    diff[2] = sum(gv27_31 - gv27_27)
    diff[3] = sum(gv27_33 - gv27_27)
    diff[4] = sum(gv27_35 - gv27_27)
    diff[5] = sum(gv27_37 - gv27_27)
    diff[6] = sum(gv27_39 - gv27_27)
    diff[7] = sum(gv27_41 - gv27_27)
    diff[8] = sum(gv27_43 - gv27_27)
    diff[9] = sum(gv27_45 - gv27_27)
    diff[10] = sum(gv27_47 - gv27_27)
    diff[11] = sum(gv27_49 - gv27_27)
    diff[12] = sum(gv27_51 - gv27_27)
    diff[13] = sum(gv27_53 - gv27_27)
    diff[14] = sum(gv27_55 - gv27_27)
    diff[15] = sum(gv27_57 - gv27_27)

# Bar chart showing the gravity changes from one timestep to the next
time1 = [27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57]
fig10, ax = plt.subplots()
plt.grid()
x = np.arange(1, len(diff) + 1)
plt.bar(x, diff, align='center')
plt.xticks(x, time1)
plt.xlabel('$Timesteps (years)$', fontsize=13)
plt.ylabel('$(\u03BCGal)$', fontsize=13)
plt.show()
fig10.savefig("overview")

# result1 = [fixed_27[0], fixed_29[0], fixed_31[0], fixed_33[0], fixed_35[0], fixed_37[0], fixed_39[0], fixed_41[0],
#            fixed_43[0], fixed_45[0], fixed_47[0],  fixed_49[0], fixed_51[0],  fixed_53[0], fixed_55[0], fixed_57[0]]
# result2 = [fixed_27[1], fixed_29[1], fixed_31[1], fixed_33[1], fixed_35[1], fixed_37[1], fixed_39[1], fixed_41[1],
#            fixed_43[1], fixed_45[1], fixed_47[1],  fixed_49[1], fixed_51[1],  fixed_53[1], fixed_55[1], fixed_57[1]]
# result3 = [fixed_27[2], fixed_29[2], fixed_31[2], fixed_33[2], fixed_35[2], fixed_37[2], fixed_39[2], fixed_41[2],
#            fixed_43[2], fixed_45[2], fixed_47[2],  fixed_49[2], fixed_51[2],  fixed_53[2], fixed_55[2], fixed_57[2]]
# result4 = [fixed_27[3], fixed_29[3], fixed_31[3], fixed_33[3], fixed_35[3], fixed_37[3], fixed_39[3], fixed_41[3],
#            fixed_43[3], fixed_45[3], fixed_47[3],  fixed_49[3], fixed_51[3],  fixed_53[3], fixed_55[3], fixed_57[3]]
# result5 = [fixed_27[4], fixed_29[4], fixed_31[4], fixed_33[4], fixed_35[4], fixed_37[4], fixed_39[4], fixed_41[4],
#            fixed_43[4], fixed_45[4], fixed_47[4],  fixed_49[4], fixed_51[4],  fixed_53[4], fixed_55[4], fixed_57[4]]


fig11, ax = plt.subplots()
x = np.arange(16)
rcParams['font.weight'] = 'bold'
# ax.plot(result1, 'k', label='Station 1')
# ax.plot(result2, 'k--', label='Station 2')
# ax.plot(result3, 'k:', label='Station 3')
# ax.plot(result4, 'k-.', label='Station 4')
# ax.plot(result5, '^k', label='Station 5')
ax.plot(result1, 'k', label='Station 1')
ax.plot(result2, 'm', label='Station 2')
ax.plot(result3, 'g', label='Station 3')
ax.plot(result4, 'r', label='Station 4')
ax.plot(result5, 'c', label='Station 5')
legend = ax.legend(loc='upper right', shadow=True, fontsize='medium')
# Put a nicer background color on the legend.
legend.get_frame()
# plt.title('Gravity changes for three measuring stations')
plt.xlabel('Timesteps (years)', fontsize=18, weight='bold')
plt.ylabel('$(\u03BC$Gal)', fontsize=18, weight='bold')
plt.xticks(x, time1)
plt.show()
fig11.savefig("final")

print('Result1', result1)
print('Result2', result2)
print('Result3', result3)
print('Result4', result4)
print('Result5', result5)
print('Diff', diff)
