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
from dsw_dpt_slices import dz, fluid_type
from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
from files_readme import sgas_27, sgas_29, sgas_31, sgas_33, sgas_35, sgas_39, sgas_45, sgas_49, sgas_53, sgas_57
from files_readme import phi, dpt

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

a = 1  # Turtuosity factor in moderately porous sand
# rw = 0.18  # Formation water resistivity
y = 2  # Cementation factor


# The function calculating and generating ATR maps
def atr_gen(sat_1, sat_2, soil_1, soil_2, sgas_1, sgas_2):
    N = len(sat_1)
    rt = np.zeros(50 * 50 * 50)

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
    for k in range(50):
        for j in range(50):
            for i in range(50):
                if sat_1[m] != 0 and sat_2[m] != 0 and phi[m] != 0:
                    if sat_2[m] < 0.1:
                        if sgas_2[m] > 0.7:
                            x = 1.8

                        else:
                            x = 2.5 - ((dpt[m] - 627) * 0.0277)
                    elif sat_2[m] == 1:
                        rt[m] = random.uniform(0.7, 1.1)
                        m += 1
                        continue
                    else:
                        x = 2.5 - ((dpt[m] - 630) * 0.0377)
                    if max_x < x:
                        max_x = x

                    atr2 = (a * rw) * (1 / (power(sat_2[m], x) * power(phi[m], y)))
                    rt[m] = atr2
                    m += 1
                else:
                    rt[m] = 0.1
                    m += 1

    Rt = np.reshape(rt, (50, 50, 50))

    result = Rt * GRID
    for k in range(50):
        for j in range(50):
            for i in range(50):
                if result[k, j, i] == 0:
                    result[k, j, i] = random.uniform(power(10, 0.8), power(10, 1.3))

    print("This is the highest value of the saturation exponent:", max_x)

    return result


########################################################################################################################


########################################################################################################################
def grid(d_z, phi):
    m = 0
    dz = np.zeros(50 * 50 * 50)
    for k in range(50):
        for j in range(50):
            for i in range(50):
                dz[m] = d_z[m]
                m += 1
    dz = np.reshape(dz, (50, 50, 50))
    poro = np.zeros((50, 50, 50))
    return dz, poro


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
GRID, POR = grid(dz, phi)

atr_2727 = atr_gen(sat27, sat27, soil_27, soil_27, sgas_27, sgas_27)
atr_2729 = atr_gen(sat27, sat29, soil_27, soil_29, sgas_27, sgas_29)
atr_2731 = atr_gen(sat27, sat31, soil_27, soil_31, sgas_27, sgas_31)
atr_2733 = atr_gen(sat27, sat33, soil_27, soil_33, sgas_27, sgas_33)
atr_2735 = atr_gen(sat27, sat35, soil_27, soil_35, sgas_27, sgas_35)
atr_2739 = atr_gen(sat27, sat39, soil_27, soil_39, sgas_27, sgas_39)
atr_2745 = atr_gen(sat27, sat45, soil_27, soil_45, sgas_27, sgas_45)
atr_2749 = atr_gen(sat27, sat49, soil_27, soil_49, sgas_27, sgas_49)
atr_2753 = atr_gen(sat27, sat53, soil_27, soil_53, sgas_27, sgas_53)
atr_2757 = atr_gen(sat27, sat57, soil_27, soil_57, sgas_27, sgas_57)

labels = [atr_2727, atr_2729, atr_2731, atr_2733, atr_2735, atr_2739, atr_2745, atr_2749, atr_2753, atr_2757]
summer = np.zeros((50, 50, 50))
data = 0
temp = 0
# for k in range(50):
#     data = labels[k]
#     for i in range(50):
#         temp += data[i]
#     summer[k] = temp / 50
#     temp = 0
########################################################################################################################

########################################################################################################################
plt.style.use('classic')

X = np.linspace(598982.08, 610735.23, 50)
Y = np.linspace(8146475.91, 8160423.65, 50)

dt_labels_tmstp = [atr_2727, atr_2729, atr_2731, atr_2733, atr_2735, atr_2739, atr_2745, atr_2749, atr_2753, atr_2757]
labels_tmstp = ["atr_2727", "atr_2729", "atr_2731", "atr_2733", "atr_2735", "atr_2739", "atr_2745", "atr_2749",
                "atr_2753", "atr_2757"]

for s in range(0, 10):
    dt = labels[s]
    outF = open("myOutFile_" + labels_tmstp[s] + ".txt", "w")
    for i in range(50):
        for j in range(50):
            outF.write(str(X[i]) + " " + str(Y[49 - j]) + " " + str(1 / dt[j, i]))
            outF.write("\n")
    outF.close()

LOGMIN = 0.1

########################################################################################################################

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

########################################################################################################################

for i in range(0, 10):
    dt = labels[i]
    for j in range(50):
        fig = plt.figure(figsize=(6.4, 4.8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        x_bronn_8_1 = 9592.92  # X coordinate of the well 7324/8-1 at cell[114, 59]
        y_bronn_8_1 = 5279  # Y coordinate of the well 7324/8-1 at cell[114, 59]
        f = RectBivariateSpline(Y * 0.001, X * 0.001, dt[j])
        Z = f(Y * 0.001, X * 0.001)
        Z = gaussian_filter(dt[j], sigma=0, mode='reflect')
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
        plt.clim(0, 4)
        plt.title(labels_tmstp[i] + "_layer " + str(j + 1) + "_" + phs)
        plt.xlabel("X [$km$]")
        plt.ylabel("Y [$km$]")
        ax = plt.gca()
        ax.patch.set_visible(False)
        fig.patch.set_visible(False)
        ax.grid(color='b', linestyle='-', linewidth=0.5)
        plt.title(labels_tmstp[i] + "_layer " + str(j + 1) + "_" + phs)
        plt.show()
        fig.savefig(labels_tmstp[i] + "__" + str(j + 1) + "_" + phs)
        ################################################################################################################
        dt = dt.ravel()
        dt_modified = np.delete(dt, np.where(dt == 0.1))
        plot_loghist(dt_modified, 37)
        mn, mx = plt.xlim()
        plt.xlim(mn, mx)
        plt.ylabel('%')
        plt.xlabel('')
        ax.set_xscale('log')
        plt.title("")
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.show()
