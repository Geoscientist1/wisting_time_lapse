import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from matplotlib import rcParams

A = 157
B = 159

# These are dense data grid dense simulation grid (101) responses (the original ones)
# grav_data_32_101_main.txt
# grav_data_42_101_main.txt
# grav_data_57_101_main.txt

# These are dense data grid (157*159) coarse simulation grid (5) responses
# grav_data_32_coarse.txt
# grav_data_42_coarse.txt
# grav_data_57_coarse.txt


# These are coarse (55*57) data grid coarse simulation grid (5) responses
# grav_data_32_coarse_dt_grid.txt
# grav_data_42_coarse_dt_grid.txt
# grav_data_57_coarse_dt_grid.txt

gv27_57 = np.loadtxt("grav_data_57_coarse.txt", delimiter=" ")
gv_real = np.loadtxt("grav_data_57_101_main.txt", delimiter=" ")
gv_dt = np.loadtxt("grav_data_32_coarse_dt_grid.txt", delimiter=" ")
gvg = np.loadtxt("fft_gv_data.txt", delimiter=" ")
gvl = np.loadtxt("delta_rho.txt", delimiter=" ")
gv_diff = np.absolute(gvl - gvg)
sat27 = np.loadtxt('sat27_coarse.txt', delimiter=' ')
sat57 = np.loadtxt('sat57_coarse.txt', delimiter=' ')
dswat_observed = np.subtract(sat57, sat27)

gt = np.loadtxt("recovered_data_temp_118_9_10000.txt", delimiter=" ")
ga = np.loadtxt("grav_data_57_coarse_dt_grid.txt", delimiter=" ")
gb = np.loadtxt("recovered_data_temp_172_3_100000.txt", delimiter=" ")
gg = np.loadtxt("grav_data_start_model_dt_grid.txt", delimiter=" ")


t1 = np.loadtxt("without_grid_expansion_FFT.txt", delimiter=" ")
t2 = np.loadtxt("with_grid_expansion_FFT.txt", delimiter=" ")

gs = np.absolute(t1 - t2)

# gv_inv = np.loadtxt("grav_data_57_101_main.txt", delimiter=" ")
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
plt.style.use('classic')

dt_labels_tmstp = [gb]
labels_tmstp = ["gv_diff"]
rcParams['font.weight'] = 'bold'
dt = dt_labels_tmstp[0]
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
# f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
# Z = f(Y * 0.001, X * 0.001)
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
# plt.clim(0, 50)
# plt.clim(0, 13)
# plt.axis('off')
plt.title("", fontsize=15, weight='bold')
plt.title("Difference", fontsize=15, weight='bold')
plt.xlabel("East [km]", fontsize=15, weight='bold')
# ax.xaxis.tick_top()
plt.ylabel("North [km]", fontsize=15, weight='bold')
ax = plt.gca()
ax.patch.set_visible(False)
fig.patch.set_visible(False)
# legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
# ax.legend()
# Put a nicer background color on the legend.
# legend.get_frame()
# ax.set_xticks(X)
# ax.set_yticks(Y)
# ax.grid(color='b', linestyle='-', linewidth=0.5)
# ax.set_xticklabels(np.arange(1, 158, 1))
# ax.set_yticklabels(np.arange(1, 160, 1))
plt.show()
# fig.savefig(labels_tmstp[i], transparent=True)
