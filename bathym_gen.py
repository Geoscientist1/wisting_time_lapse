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
from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
from files_readme import phi, dpt


########################################################################################################################

# This script models and generates a bathymetry surface to be used in the forward modeling in sblwiz
# The script also plots the modeled bathymetry  surface

########################################################################################################################

# The function calculating and generating ATR maps
def bath_gen():
    rt = np.zeros(k_range * 157 * 159)
    m = 0
    for k in range(k_range):
        for j in range(159):
            for i in range(157):
                rt[m] = random.uniform(power(10, 0.1), power(10, 0.5))
                m += 1

    Rt = np.reshape(rt, (k_range, 159, 157))

    if typo == "stk":
        stk_atr = np.zeros((159, 157))
        temp = 0
        for i in range(157):
            for j in range(159):
                for k in range(k_range):
                    temp += (Rt[k, j, i])
                stk_atr[j, i] = temp
                temp = 0
        result = stk_atr
        for j in range(159):
            for i in range(157):
                # if result[j, i] < 0.1:
                if result[j, i] == 0:
                    result[j, i] = random.uniform(power(10, 0.1), power(10, 0.5))
                # result[j, i] = random.uniform(power(10, 0.8), power(10, 1.3))
                # if result[j, i] == 0.1:
                #     result[j, i] = bck_res_compress[j, i]

    elif typo == "NONE":
        result = Rt
        for k in range(k_range):
            for j in range(159):
                for i in range(157):
                    result[k, j, i] = random.uniform(power(10, 0.7), power(10, 0.7))

    dRt = np.zeros((10, 159, 157))

    return result


########################################################################################################################
bath_srfc = bath_gen()
########################################################################################################################
plt.style.use('classic')

X = np.linspace(598982.08, 610735.23, 157)
Y = np.linspace(8146475.91, 8160423.65, 159)
Z = np.linspace(-1082.86, -596.02, k_range)
########################################################################################################################
# This part of the script is writing out the modeled data into a Gslib friendly format to be used in SBLwiz software
outF = open("myOutFile_" + "bath_srfc" + ".txt", "w")
for k in range(k_range):
    for j in range(159):
        for i in range(157):
            if bath_srfc[k, j, i] <= 0:
                outF.write("{0} {1} {2} {3} {4} {5} {6}".format(str(i + 1), str(159 - j), str(k + 1),
                                                                str(X[i]), str(Y[j]), str(Z[k]),
                                                                str(-99999.00)))
            else:
                outF.write("{0} {1} {2} {3} {4} {5} {6}".format(str(i + 1), str(159 - j), str(k + 1),
                                                                str(X[i]), str(Y[j]), str(Z[k]),
                                                                str(bath_srfc[k, j, i])))
            outF.write("\n")
outF.close()

LOGMIN = 0.1

########################################################################################################################


if typo == "stk":
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')
    f = RectBivariateSpline(Y * 0.001, X * 0.001, bath_srfc)
    x_bronn_8_1 = 9592.92  # X coordinate of the well 7324/8-1 at cell[114, 59]
    y_bronn_8_1 = 5279  # Y coordinate of the well 7324/8-1 at cell[114, 59]
    # plt.plot([x_bronn_8_1], [y_bronn_8_1], marker='o', markersize=7, color="k")
    Z = f(Y * 0.001, X * 0.001)
    Z = gaussian_filter(bath_srfc, sigma=0, mode='reflect')
    plt.imshow(np.log10(Z), extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
               cmap=plt.cm.get_cmap('jet', 1000),
               interpolation='nearest', origin='upper')

    clb = plt.colorbar()
    tick_locator = ticker.MaxNLocator(nbins=7)
    clb.locator = tick_locator
    clb.update_ticks()
    clb.set_label(r'log10(ATR)[$\Omega m^{2}$]', rotation=90, fontsize=15)
    ax.set_ylabel("Y(m)", labelpad=15)
    plt.clim(0, 4)
    # plt.title(labels_tmstp[i] + "_layerstack_" + phs)
    plt.xlabel("X [$km$]")
    # ax.xaxis.tick_top()
    plt.ylabel("Y [$km$]")
    ax = plt.gca()
    ax.patch.set_visible(False)
    fig.patch.set_visible(False)
    ax.grid(color='b', linestyle='-', linewidth=0.5)
    plt.title("Bathymetry Surface")
    plt.show()
    fig.savefig("bath_srfc")  # + ".pdf")
########################################################################################################################
else:
    for j in range(k_range):
        fig = plt.figure(figsize=(6.4, 4.8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        f = RectBivariateSpline(Y, X, bath_srfc[j])
        Z = f(Y, X)
        Z = gaussian_filter(bath_srfc[j], sigma=0.27)
        plt.imshow(np.log10(Z), extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
                   interpolation='spline36', origin='upper')
        clb = plt.colorbar()
        clb.set_label(r'log10(ATR)[$\Omega m^{2}$]', rotation=90, fontsize=15)
        ax.set_ylabel("Y(m)", labelpad=15)
        plt.clim(0, 5)
        plt.title("bath_srfc" + "_layer " + str(j + 1))
        plt.xlabel("X [$m$]")
        ax.xaxis.tick_top()
        plt.ylabel("Y [$m$]")
        plt.show()
        fig.savefig("bath_srfc" + "__" + str(j + 1))
