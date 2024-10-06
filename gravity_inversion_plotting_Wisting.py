import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import scipy.interpolate as si
import pandas as pd
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
import time
# from Grav_Modeling_coarse_grid import *
from gravity_inversion_Wisting import *

X = np.linspace(0, x_max, x_range_data)
Y = np.linspace(0, y_max, y_range_data)
plt.style.use('classic')

########################################################################################################################
dt = np.loadtxt("recovered_data_625.txt", delimiter=" ")
labels_tmstp = ["Initial data"]
########################################################################################################################
# dt = dt_labels_tmstp[i]
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
clb.ax.set_title('$(\u03BC$Gal)', fontsize=15)
ax.set_ylabel("Y(m)", labelpad=15)
# if labels_tmstp[i] == "recovered_datad":
#     print("")
#     # plt.clim(0.04, 0.08)
# elif labels_tmstp[i] == "difference":
#     print("")
#     # plt.clim(0.04, 0.08)
# else:
#     plt.clim(0.05, 0.27)
# plt.axis('off')
# plt.title(labels_tmstp[i] + "_layerstack_" + phs)
# plt.title(str(labels_tmstp[i]))
plt.xlabel("X [$km$]")
# ax.xaxis.tick_top()
plt.ylabel("Y [$km$]")
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
# plt.clim(0, 20)
plt.show()
# fig.savefig(labels_tmstp[i], transparent=True)
