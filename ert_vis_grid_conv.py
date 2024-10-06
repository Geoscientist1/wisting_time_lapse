import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from scipy import interpolate
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import copy
from matplotlib import colors
from matplotlib import ticker
from matplotlib.patches import Rectangle
from numpy import *
import matplotlib.cm as cm
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline


def ert_vis_grid_conv(resis_map, cell_hgt):
    line_count = 0
    resist = np.zeros((318, 314))
    height = np.zeros((318, 314))
    wat_sat = np.zeros((318, 314))
    bck_res = np.zeros((318, 314))

    with open(resis_map, "r") as file:
        for line in file:
            if line != "\n":
                line_count += 1

            arr = line.split()
            for k in range(len(arr)):
                i = int(arr[3])
                j = 318 - int(arr[4])
                resist[j, i] = float(arr[2])

    with open(cell_hgt, "r") as file:
        for line in file:
            if line != "\n":
                line_count += 1

            arr1 = line.split()
            for k in range(len(arr1)):
                i = int(arr1[3])
                j = 318 - int(arr1[4])
                height[j, i] = float(arr1[2])

    with open(wat_sat_27, "r") as file:
        for line in file:
            if line != "\n":
                line_count += 1

            arr2 = line.split()
            for k in range(len(arr2)):
                i = int(arr2[3])
                j = 318 - int(arr2[4])
                wat_sat[j, i] = float(arr2[2])

    resist = resist * height

    for j in range(318):
        for i in range(314):
            if resist[j, i] == 0:
                # or wat_sat[j, i] == 1:
                resist[j, i] = 0.1
                # resist[j, i] = random.uniform(power(10, 0.8), power(10, 1.3))
            if wat_sat[j, i] == 1:
                bck_res[j, i] = resist[j, i]

    return resist, bck_res


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
    labels = ["ATR_csem_inversion"]
    plt.legend(handles, labels)


########################################################################################################################


########################################################################################################################
atr_csem = 'atr_csem_inversion.txt'
atr_mod_27 = 'atr_mod_27_t.txt'
wat_sat_27 = 'wat_sat_27.txt'

dz = 'dzmap.txt'

resist, bck_res = ert_vis_grid_conv(atr_csem, dz)
resist_mod_27, bck_res_mod_27 = ert_vis_grid_conv(atr_mod_27, dz)

plt.style.use('classic')

X = np.linspace(598982.08, 610735.23, 314)
Y = np.linspace(8146475.91, 8160423.65, 318)

X = np.linspace(0, 11367.59, 314)
Y = np.linspace(0, 11557.04, 318)


outF = open("myOutFile_" + "atr_csem_inversion" + ".txt", "w")
for i in range(314):
    for j in range(318):
        outF.write(str(X[i]) + " " + str(Y[317 - j]) + " " + str(1 / resist[j, i]))
        outF.write("\n")
outF.close()

LOGMIN = 0.1
########################################################################################################################
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
f = RectBivariateSpline(Y * 0.001, X * 0.001, resist)
x_bronn_8_1 = 9.46  # X coordinate of the well 7324/8-1 at cell[114, 59]
y_bronn_8_1 = 7.59  # Y coordinate of the well 7324/8-1 at cell[114, 59]
plt.plot([x_bronn_8_1], [y_bronn_8_1], marker='o', markersize=7, color="w", linestyle='None', label="Well 7324/8-1")
Z = f(Y * 0.001, X * 0.001)
Z = gaussian_filter(resist, sigma=0, mode='reflect')
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
# plt.title(labels_tmstp[i] + "_layerstack_" + phs)
plt.xlabel("X [$km$]")
# ax.xaxis.tick_top()
plt.ylabel("Y [$km$]")
ax = plt.gca()
# ax.set_xticks(X)
# ax.set_yticks(Y)
ax.grid(color='b', linestyle='-', linewidth=0.5)
legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
legend.get_frame()
# ax.set_xticklabels(np.arange(1, 158, 1))
# ax.set_yticklabels(np.arange(1, 160, 1))
plt.show()
fig.savefig("atr_csem_data")
########################################################################################################################
resist = resist.ravel()
resist_modified = np.delete(resist, np.where(resist == 0.1))
plot_loghist(resist_modified, 37)
mn, mx = plt.xlim()
plt.xlim(mn, mx)
plt.ylabel('%')
plt.xlabel('')
ax.set_xscale('log')
plt.title("")
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.show()
########################################################################################################################

########################################################################################################################


fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
f = RectBivariateSpline(Y, X, resist_mod_27)
x_bronn_8_1 = 9592.92 + 598982.08  # X coordinate of the well 7324/8-1 at cell[114, 59]
y_bronn_8_1 = 5279 + 8146475.91  # Y coordinate of the well 7324/8-1 at cell[114, 59]
plt.plot([x_bronn_8_1], [y_bronn_8_1], marker='o', markersize=7, color="k")
Z = f(Y, X)
Z = gaussian_filter(resist_mod_27, sigma=0, mode='reflect')
plt.imshow(np.log10(Z), extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
           # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
           interpolation='nearest', origin='upper')

clb = plt.colorbar()
clb.set_label(r'log10(ATR)[$\Omega m^{2}$]', rotation=90, fontsize=15)
ax.set_ylabel("Y(m)", labelpad=15)
plt.clim(0, 4)
# plt.title(labels_tmstp[i] + "_layerstack_" + phs)
plt.xlabel("X [$m$]")
ax.xaxis.tick_top()
plt.ylabel("Y [$m$]")
ax = plt.gca()
# ax.set_xticks(X)
# ax.set_yticks(Y)
ax.grid(color='b', linestyle='-', linewidth=0.5)
# ax.set_xticklabels(np.arange(1, 158, 1))
# ax.set_yticklabels(np.arange(1, 160, 1))
plt.show()
# fig.savefig(labels_tmstp[i] + "_layerstack_" + phs + ".pdf")
########################################################################################################################
resist_mod_27 = resist_mod_27.ravel()
resist_mod_27_modified = np.delete(resist_mod_27, np.where(resist == 0.1))
plot_loghist(resist_mod_27_modified, 37)
mn, mx = plt.xlim()
plt.xlim(mn, mx)
plt.ylabel('%')
plt.xlabel('')
ax.set_xscale('log')
plt.title("")
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.show()
