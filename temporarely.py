import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from scipy import interpolate
from scipy.ndimage import zoom
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


def ert_vis_grid_conv(resis_map, cell_hgt):
    line_count = 0
    resist = np.zeros((318, 314))
    height = np.zeros((318, 314))

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

    resist = resist * height

    for j in range(318):
        for i in range(314):
            if resist[j, i] == 0:
                resist[j, i] = random.uniform(power(10, 0.8), power(10, 1.3))

    return resist


atr_csem = 'atr_csem_inversion.txt'
dz = 'dzmap.txt'
resist = ert_vis_grid_conv(atr_csem, dz)
resist = zoom(resist, (0.17, 0.17))

for j in range(54):
    for i in range(53):
        if resist[j, i] < 0:
            resist[j, i] = 1000

plt.style.use('classic')

X = np.linspace(598982.08, 610735.23, 53)
Y = np.linspace(8146475.91, 8160423.65, 54)
Z = np.linspace(0, 278, 50)

outF = open("myOutFile_" + "atr_csem_inversion" + ".txt", "w")
temp = 0
for k in range(50):
    for i in range(53):
        for j in range(54):
            # outF.write(str(X[i]) + " " + str(Y[52 - j]) + " " + str(1 / resist[j, i]))
            outF.write(str(X[i]) + " " + str(Y[52 - j]) + " " + str(temp) + " " + str(1 / resist[j, i]))
            outF.write("\n")
    temp += 5
outF.close()

LOGMIN = 0.1
########################################################################################################################
# Just to know the number of lines we have in the file
with open("myOutFile_atr_csem_inversion.txt", "r") as f:
    print(len(f.readlines()))

########################################################################################################################
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
f = RectBivariateSpline(Y * 0.001, X * 0.001, resist)
x_bronn_8_1 = 10.9326  # X coordinate of the well 7324/8-1 at cell[114, 59]
y_bronn_8_1 = 9.73446  # Y coordinate of the well 7324/8-1 at cell[114, 59]
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
