import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from scipy import interpolate
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import copy
from matplotlib import colors
from matplotlib.patches import Rectangle
from matplotlib import ticker
from numpy import *
import matplotlib.cm as cm
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
import scipy.stats as st
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from dsw_dpt_slices import GRID, POR, dz, k_range, fluid_type, typo
from files_readme import sat27, sat29, sat31, sat32, sat33, sat35, sat37, sat39, sat42, sat45, sat47, sat49, sat52, \
    sat53, sat57
from files_readme import soil_27, soil_29, soil_31, soil_32, soil_33, soil_35, soil_37, soil_39, soil_42, soil_45, \
    soil_47, soil_49, soil_52, soil_53, soil_57
from files_readme import sgas_27, sgas_29, sgas_31, sgas_32, sgas_33, sgas_35, sgas_37, sgas_39, sgas_42, sgas_45, \
    sgas_47, sgas_49, sgas_52, sgas_53, sgas_57
from files_readme import phi, dpt
from files_readme import A, B

########################################################################################################################

# This script models and generates ATR maps from Petrel simulation grids.
# It calculates the formation resistivity Rt at each cell, and then
# integrates the product of cell resistivity and cell height across the
# whole reservoir.
# The script also plots the resulted ATR maps accompanied with their
# corresponding histograms.

########################################################################################################################
rho_oil = 835  # Oil density taken from NPD factpages
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
typo = "stk"

a = 1  # Turtuosity factor in moderately porous sand
# rw = 0.18  # Formation water resistivity
y = 2  # Cementation factor


# The function calculating and generating ATR maps
def atr_gen(sat_1, sat_2, soil_1, soil_2, sgas_1, sgas_2):
    rt = np.zeros(k_range * A * B)

    if fluid_type == "Water-Oil_saturation":
        sat_1 = 1 - soil_1
        sat_2 = 1 - soil_2
    elif fluid_type == "Water-Gas_saturation":
        sat_1 = 1 - sgas_1
        sat_2 = 1 - sgas_2
    elif fluid_type == "Water_saturation":
        sat_1 = sat_1
        sat_2 = sat_2

    rw = 0.18
    m = 0
    max_x = 0
    for k in range(k_range):
        for j in range(B):
            for i in range(A):
                if sat_1[m] != 0 and sat_2[m] != 0 and phi[m] != 0:

                    # if sgas_2[m] > 0.3:
                    #     x = 1.8
                    #     if sat_2[m] < 0.1:
                    #         rw = 1
                    #         atr2 = (a * rw) * (1 / (power(sat_2[m], x) * power(phi[m], y)))
                    #         rt[m] = atr2
                    #         m += 1
                    #         continue
                    #     else:
                    #         rw = 0.18
                    #         atr2 = (a * rw) * (1 / (power(sat_2[m], x) * power(phi[m], y)))
                    #         rt[m] = atr2
                    #         m += 1
                    #         continue
                    # if sat_2[m] < 0.1:
                    #     if sat_2[m] < 0.08:
                    #         x = 3.15
                    #         rw = 1
                    #     else:
                    #         x = (-79.69 * power(sat_2[m], 3)) + (56.76 * power(sat_2[m], 2)) + (
                    #             -14 * power(sat_2[m], 1)) + 3.59
                    #     # x = 3.15
                    #         rw = 1
                    # elif sgas_2[m] < 0.3 and sat_2[m] > 0.1:
                    #     # x = (-79.69 * power(sat_2[m], 3)) + (56.76 * power(sat_2[m], 2)) + (
                    #     #         -14 * power(sat_2[m], 1)) + 3.59
                    #     x = 2.1
                    #     rw = 0.18
                    if sat_2[m] < 0.1:
                        if sgas_2[m] > 0.7:
                            x = 1.8

                        else:
                            x = 2.5 - ((dpt[m] - 627) * 0.0277)
                            # if x < 0:
                            #     x = 0.1
                    elif sat_2[m] == 1:
                        # rt[m] = random.uniform(0.7, 1.1)
                        rt[m] = 0.13
                        # rt[m] = bck_res[m]
                        # rt[m] = 0
                        m += 1
                        continue
                    else:
                        x = 2.5 - ((dpt[m] - 630) * 0.0377)
                        # if x < 0:
                        #     x = 0.1
                        # if sat_2[m] == 1:
                        #     rw_temp = 0.1
                        #     rt[m] = (a * rw_temp) * (1 / (power(sat_2[m], x) * power(phi[m], y)))
                        #     m += 1
                        #     continue

                    # # atr1 = (a * rw) * (1 / (power(sat_1[m] / 10, x) * power(phi[m], y)))
                    # print("This is x:", x)
                    if max_x < x:
                        max_x = x

                    atr2 = (a * rw) * (1 / (power((sat_2[m]), x) * power(phi[m], y)))
                    atr1 = (a * rw) * (1 / (power((sat_1[m]), x) * power(phi[m], y)))
                    atr = np.abs(atr1 - atr2) / atr1
                    rt[m] = atr  # 7 ohm.m is background resistivity
                    if rt[m] <= 0:
                        rt[m] = 0.1
                    m += 1
                else:
                    rt[m] = 0.1
                    m += 1

    Rt = np.reshape(rt, (k_range, B, A))

    if typo == "stk":
        stk_atr = np.zeros((B, A))
        # bck_res_ref = bck_res.flatten()
        # bck_res_interp = interpolate.interp1d(np.arange(bck_res.size), bck_res_ref)
        # bck_res_compress = bck_res_interp(np.linspace(0, bck_res.size - 1, stk_atr.size))
        temp = 0
        for i in range(A):
            for j in range(B):
                for k in range(k_range):
                    temp += (Rt[k, j, i] * GRID[k, j, i])
                stk_atr[j, i] = temp
                temp = 0
        result = stk_atr
        for j in range(B):
            for i in range(A):
                if result[j, i] < 0.1:
                    # if result[j, i] == 0:
                    result[j, i] = 0.1
                # result[j, i] = random.uniform(power(10, 0.8), power(10, 1.3))
                # if result[j, i] == 0.1:
                #     result[j, i] = bck_res_compress[j, i]

    elif typo == "NONE":
        result = Rt  # the input model in sblwiz is a resistivity model, so we do not need to multiply the
        # resistivity with the grid thickness
        # result = Rt * GRID
        for k in range(k_range):
            for j in range(B):
                for i in range(A):
                    if result[k, j, i] <= 0.1:
                        result[k, j, i] = random.uniform(power(10, 0.9), power(10, 1.2))
                        # result[k, j, i] = random.uniform(power(10, 0.8), power(10, 1.3))
                    # elif result[k, j, i] == 0:
                    #     result[k, j, i] = random.uniform(power(10, 0.8), power(10, 1.3))

    dRt = np.zeros((10, B, A))
    # s = 0
    # fk = 1
    # temp = 0
    # for i in range(A):
    #     for j in range(B):
    #         for n in range(101):
    #             temp += Rt[n, j, i]
    #             if n != 0 and mod(n, 9) == 0 and n != 99:
    #                 dRt[s, j, i] = temp / (fk * 10)
    #                 s += 1
    #                 temp = 0  # inside the loop or outside depends on whether we need the accumulated average or
    #         # single layer average
    #         s = 0
    #         fk = 1

    return result


########################################################################################################################


########################################################################################################################
# The function calculating water saturation along a well trajectory
def Sw_profile(sat_1, dz):
    Sw = np.zeros(101)
    tvd = np.zeros(101)
    cell_x = 114
    cell_y = 59
    m = 0
    n = 0
    for k in range(101):
        for j in range(B):
            for i in range(A):
                if i == cell_x and j == cell_y and dz[m] != 0:
                    Sw[n] = sat_1[m]
                    tvd[n] = dz[m]
                    n += 1
                m += 1
    tvd = np.cumsum(tvd)
    return Sw, tvd


########################################################################################################################
# The function plotting the histograms
def plot_loghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    a = bins[0]
    b = bins[-1]
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    N, bins, patches = plt.hist(x, bins=logbins, color='tab:blue', weights=np.ones(len(x)) / len(x))
    plt.xscale('log')
    # create legend
    cmap = plt.get_cmap('tab10')
    Color = cmap(0)
    for i in range(0, 4):
        patches[i].set_facecolor(Color)
    handles = [Rectangle((0, 0), 1, 1, color=Color, ec="k")]
    labels = ["ATR_mod_27_t"]
    plt.legend(handles, labels)


########################################################################################################################
atr_2727 = atr_gen(sat27, sat27, soil_27, soil_27, sgas_27, sgas_27)
atr_2732 = atr_gen(sat27, sat32, soil_27, soil_32, sgas_27, sgas_32)
atr_2737 = atr_gen(sat27, sat37, soil_27, soil_37, sgas_27, sgas_37)
atr_2742 = atr_gen(sat27, sat42, soil_27, soil_42, sgas_27, sgas_42)
atr_2747 = atr_gen(sat27, sat47, soil_27, soil_47, sgas_27, sgas_47)
atr_2752 = atr_gen(sat27, sat52, soil_27, soil_52, sgas_27, sgas_52)
atr_2757 = atr_gen(sat27, sat57, soil_27, soil_57, sgas_27, sgas_57)

labels = [atr_2727, atr_2732, atr_2737, atr_2742, atr_2747, atr_2752, atr_2757]
summer = np.zeros((10, B, A))
data = 0
temp = 0
for k in range(len(labels)):
    data = labels[k]
    for i in range(k_range):
        temp += data[i]
    summer[k] = temp / 10
    temp = 0
########################################################################################################################
Sw_27, TVD = Sw_profile(sat27, dz)
# Sw_29 = Sw_profile(sat29, dz)
# Sw_31 = Sw_profile(sat29, dz)
# Sw_27 = Sw_profile(sat27, dz)
# Sw_27 = Sw_profile(sat27, dz)
# Sw_27 = Sw_profile(sat27, dz)
# Sw_27 = Sw_profile(sat27, dz)
# Sw_27 = Sw_profile(sat27, dz)
# Sw_27 = Sw_profile(sat27, dz)
# Sw_27 = Sw_profile(sat27, dz)

# fig = plt.figure(facecolor='#e1ddbf')
# ax = fig.add_subplot(111)
# plt.plot(Sw_27, TVD, label='Sw', color='c', linewidth=3, marker='h', markerfacecolor='lightgreen',
#          markeredgewidth=1,
#          markersize=7, markevery=3)
# plt.title('Water saturation along Well 7324/8-1')
# plt.xlabel("Sw (Water Saturation)", fontsize=20)
# plt.ylabel("TVD [m]", fontsize=20)
# plt.gca().invert_yaxis()
# plt.legend(loc='best')
# plt.setp(ax.spines.values(), linewidth=3)
# # The ticks
# ax.xaxis.set_tick_params(width=3)
# ax.yaxis.set_tick_params(width=3)
# plt.locator_params(axis='x', nbins=3)
# plt.show()
# fig.savefig('Water_saturation_Profile')

########################################################################################################################
plt.style.use('classic')

X = np.linspace(0, 11367.59, A)
Y = np.linspace(0, 11557.04, B)

# X = np.linspace(599353.72, 606792.17, A)
# Y = np.linspace(8146346.56, 8160662.48, B)
Z = np.linspace(-1082.86, -596.02, k_range)

dt_labels_tmstp = [atr_2727, atr_2732, atr_2737, atr_2742, atr_2747, atr_2752, atr_2757]
labels_tmstp = ["atr_2727", "atr_2732", "atr_2737", "atr_2742", "atr_2747", "atr_2752", "atr_2757"]

# for s in range(0, 10):
#     dt = labels[s]
#     outF = open("myOutFile_" + labels_tmstp[s] + ".txt", "w")
#     for i in range(A):
#         for j in range(B):
#             outF.write(str(X[i]) + " " + str(Y[B-1 - j]) + " " + str(1 / dt[j, i]))
#             outF.write("\n")
#     outF.close()

########################################################################################################################
# This part of the script is writing out the modeled data into a Gslib friendly format to be used in SBLwiz software
for s in range(0, len(dt_labels_tmstp)):
    dt = labels[s]
    outF = open("myOutFile_" + labels_tmstp[s] + ".xyz", "w")
    # for k in range(k_range):
    for j in range(B):
        for i in range(A):
            if dt[j, i] == 0.1:
                # This parcel of the code can be used when wants the output file to be in gslib format.
                # outF.write("{0} {1} {2} {3} {4} {5} {6}".format(str(i + 1), str(B - j), str(k + 1),
                #                                                 str(X[i]), str(Y[j]), str(Z[k]),
                #                                                 str(random.uniform(power(10, 0.9),
                #                                                                    power(10, 1.2)))))
                # outF.write("{0} {1} {2}".format(str(X[i]), str(Y[j]), str(dt[j, i])))
                continue
            else:
                # outF.write("{0} {1} {2} {3} {4} {5} {6}".format(str(i + 1), str(B - j), str(k + 1),
                #                                                 str(X[i]), str(Y[j]), str(Z[k]), str(dt[j, i])))
                # This parcel of the code can be used when wants the output file to be in gslib format.
                outF.write("{0} {1} {2}".format(str(X[i]), str(Y[j]), str(dt[j, i])))
            outF.write("\n")
    outF.close()

LOGMIN = 0.1

########################################################################################################################
# for i in range(0, 10):
#     dt = dt_labels_tmstp[i]
#     for j in range(10):
#         fig = plt.figure(figsize=(6.4, 4.8))
#         ax = fig.add_subplot(111)
#         ax.set_aspect('equal')
#         f = RectBivariateSpline(Y, X, dt[j] / 10)
#         Z = f(Y, X)
#         Z = gaussian_filter(dt[j], sigma=0.4)
#         plt.imshow(np.log10(Z), extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
#                    # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
#                    interpolation='spline36', origin='upper')
#         clb = plt.colorbar()
#         clb.set_label(r'log10(ATR)[$\Omega m^{2}$]', rotation=90, fontsize=15)
#         ax.set_ylabel("Y(m)", labelpad=15)
#         k = [-0.15, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0]
#         plt.clim(-2, 7)
#         plt.title('atr ---> ' + fig_labels_tmstp[i] + ' ---> ' + fig_labels_layer[j])
#         plt.xlabel("X[m]")
#         ax.xaxis.tick_top()
#         plt.ylabel("Y[m]")
#         plt.show()
#         fig.savefig(oupt_labels_tmstp[i] + oupt_labels_layer[j])
########################################################################################################################
if fluid_type == "Water_saturation":
    phs = "w"
elif fluid_type == "Oil_saturation":
    phs = "o"
elif fluid_type == "Gas_saturation":
    phs = "g"
elif fluid_type == "Water-Oil_saturation":
    phs = "wo"
elif fluid_type == "Water-Gas_saturation":
    phs = "wg"

if typo == "stk":
    for i in range(0, len(dt_labels_tmstp)):
        dt = labels[i]
        fig = plt.figure(figsize=(6.4, 4.8))
        ax = fig.add_subplot(111)
        ax.set_aspect('auto')
        f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
        x_bronn_8_1 = 9592.92  # X coordinate of the well 7324/8-1 at cell[114, 59]
        y_bronn_8_1 = 5279  # Y coordinate of the well 7324/8-1 at cell[114, 59]
        # plt.plot([x_bronn_8_1], [y_bronn_8_1], marker='o', markersize=7, color="k")
        Z = f(Y * 0.001, X * 0.001)
        Z = gaussian_filter(dt, sigma=0, mode='reflect')
        plt.imshow(np.log10(Z), extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
                   cmap=plt.cm.get_cmap('jet', 1000),
                   # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
                   interpolation='nearest', origin='upper')

        clb = plt.colorbar()
        tick_locator = ticker.MaxNLocator(nbins=7)
        clb.locator = tick_locator
        clb.update_ticks()
        clb.set_label(r'log10(ATR)[$\Omega m^{2}$]', rotation=90, fontsize=15)
        ax.set_ylabel("Y(m)", labelpad=15)
        plt.clim(0, 1)
        # plt.title(labels_tmstp[i] + "_layerstack_" + phs)
        plt.xlabel("X [$km$]")
        # ax.xaxis.tick_top()
        plt.ylabel("Y [$km$]")
        ax = plt.gca()
        ax.patch.set_visible(False)
        fig.patch.set_visible(False)
        # ax.set_xticks(X)
        # ax.set_yticks(Y)
        ax.grid(color='b', linestyle='-', linewidth=0.5)
        # ax.set_xticklabels(np.arange(1, 158, 1))
        # ax.set_yticklabels(np.arange(1, 160, 1))
        plt.title(labels_tmstp[i])
        plt.show()
        fig.savefig(labels_tmstp[i] + "_layerstack_" + phs)  # + ".pdf")
        ################################################################################################################
        # dt = dt.ravel()
        # dt_modified = np.delete(dt, np.where(dt == 0.1))
        # plot_loghist(dt_modified, 37)
        # mn, mx = plt.xlim()
        # plt.xlim(mn, mx)
        # plt.ylabel('%')
        # plt.xlabel('')
        # ax.set_xscale('log')
        # plt.title("")
        # plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        # plt.show()
else:
    for i in range(0, len(dt_labels_tmstp)):
        dt = labels[i]
        for j in range(k_range):
            fig = plt.figure(figsize=(6.4, 4.8))
            ax = fig.add_subplot(111)
            ax.set_aspect('equal')
            f = RectBivariateSpline(Y, X, dt[j])
            Z = f(Y, X)
            Z = gaussian_filter(dt[j], sigma=0.27)
            plt.imshow(np.log10(Z), extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
                       # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
                       interpolation='spline36', origin='upper')
            clb = plt.colorbar()
            clb.set_label(r'log10(ATR)[$\Omega m^{2}$]', rotation=90, fontsize=15)
            ax.set_ylabel("Y(m)", labelpad=15)
            plt.clim(0, 5)
            plt.title(labels_tmstp[i] + "_layer " + str(j + 1) + "_" + phs)
            plt.xlabel("X [$m$]")
            ax.xaxis.tick_toptop()
            plt.ylabel("Y [$m$]")
            plt.show()
            fig.savefig(labels_tmstp[i] + "__" + str(j + 1) + "_" + phs)
