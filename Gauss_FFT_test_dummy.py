import math
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import time
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from matplotlib import rcParams

start_time = time.time()

########################################################################################################################
M = 7  # This is the # of grid cells in the x-direction
m = np.arange(0, M, 1)
delta_x = (13213.48 - 0) / M
x0 = 0
xm = x0 + (m * delta_x)
delta_kx = ((2 * np.pi) / (M * delta_x))
kxp = np.fft.fftfreq(M, delta_x) * 2 * np.pi
p = np.arange(-M / 2, (M / 2), 1)
kxp1 = p * delta_kx
IX = 157
# kx = (2 * np.pi) / M
########################################################################################################################
N = 7  # This is the # of grid cells in the y-direction
n = np.arange(0, N, 1)
delta_y = (14173.40 - 0) / N
y0 = 0
yn = y0 + (n * delta_y)
delta_ky = ((2 * np.pi) / (N * delta_y)) * 2 * np.pi
kyq = np.fft.fftfreq(N, delta_y)
IY = 159
# ky = (2 * np.pi) / N
########################################################################################################################
L = 7  # This is the # of grid cells in the z-direction
l = np.arange(0, L, 1)
delta_z = (-1082.86 - (-596.02)) / L
z0 = -596.02
zl = z0 + (l * delta_z)
delta_kz = ((2 * np.pi) / (L * delta_z))
kzw = np.fft.fftfreq(L, delta_z) * 2 * np.pi
IZ = 5
# kz = (2 * np.pi) / L
########################################################################################################################
# Gaussian quadrature nodes and weights
n = 4
x, wi = np.polynomial.legendre.leggauss(n)
nodes = np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3, -1).T
wei = (wi * wi[:, None]).ravel()
weights = (wei * wi[:, None]).ravel()
########################################################################################################################
G = 6.67 * 1e-11  # Gravity constant
phi = 0.3
XSCALE = 188
YSCALE = 202
dx_accumulated = np.linspace(0, 1321, 7)
dy_accumulated = np.linspace(0, 1417, 7)
cell_dtp = np.random.uniform(low=600, high=1000, size=(L*N*M,))
dx = np.random.uniform(low=1, high=4, size=(L*N*M,))
dy = np.random.uniform(low=1, high=4, size=(L*N*M,))
dz = np.random.uniform(low=1, high=2, size=(L*N*M,))
########################################################################################################################
delta_g = np.zeros((7, 7))
delta_rhos = np.random.uniform(low=0, high=0.7, size=(L*N*M,))
a = G * (cell_dtp - 400) * delta_rhos * (1040 - 835) * phi
b = dx * dy * dz  # this corresponds to V v(volume) in the formula for calculating delta gz
c = a * b
for y_step in range(0, 7):
    for x_step in range(0, 7):
        m = 0
        print("y", y_step)
        print("x", x_step)
        for k in range(7):
            for j in range(7):
                for i in range(7):
                    r = np.sqrt(
                        ((x_step * XSCALE) - dx_accumulated[i]) ** 2 + ((y_step * YSCALE) - dy_accumulated[j]) ** 2)
                    d = pow((pow(r, 2) + pow(cell_dtp[m] - 400, 2)), 1.5)
                    delta_g[y_step, x_step] += ((c[m]) / d)
                    m += 1
delta_g = delta_g * 1e8
########################################################################################################################
delta_rho = np.reshape(delta_rhos, (L, N, M))
fft_delta_rho = np.fft.fftn(delta_rho, axes=(1, 2))
########################################################################################################################
k = np.sqrt((np.power(kyq, 2)) + (np.power(kxp, 2)))
f1 = (1j * 4 * np.pi * G) / (np.power(k, 2) * delta_x * delta_y)
f2 = (2 * (np.sin(0.5 * kxp * delta_x)) / kxp) * (2 * (np.sin(0.5 * kyq * delta_y)) / kyq)
f3 = np.multiply(f1, f2)
f3[0] = 0.000001
result = f3
delta_g_hat = np.multiply(result, fft_delta_rho)
delta_g_fft = np.fft.ifftn(delta_g_hat, axes=(1, 2))
delta_g_2 = np.sum(delta_g_fft, axis=2)
slice_index = 0
# np.savetxt('fft_gv_data.txt', delta_g.real[2], delimiter=' ', header='Gravity data Coarse Grid for year 2057')
########################################################################################################################
########################################################################################################################
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
plt.style.use('classic')
dt_labels_tmstp = np.real(fft_delta_rho)[slice_index]
labels_tmstp = ["gv_diff"]
rcParams['font.weight'] = 'bold'
dt = dt_labels_tmstp
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
Z = gaussian_filter(dt, sigma=0, mode='reflect')
plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
           cmap=plt.cm.get_cmap('jet', 1000),
           # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
           interpolation='nearest', origin='upper')
clb = plt.colorbar()
tick_locator = ticker.MaxNLocator(nbins=7)
clb.locator = tick_locator
clb.update_ticks()
plt.title("", fontsize=15, weight='bold')
plt.title("rho_hat", fontsize=15, weight='bold')
plt.xlabel("East [km]", fontsize=15, weight='bold')
plt.ylabel("North [km]", fontsize=15, weight='bold')
ax = plt.gca()
ax.patch.set_visible(False)
fig.patch.set_visible(False)
plt.show()
########################################################################################################################
########################################################################################################################
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
plt.style.use('classic')
dt_labels_tmstp = delta_g
rcParams['font.weight'] = 'bold'
dt = dt_labels_tmstp
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
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
plt.title("", fontsize=15, weight='bold')
plt.title("Observed data", fontsize=15, weight='bold')
plt.xlabel("East [km]", fontsize=15, weight='bold')
plt.ylabel("North [km]", fontsize=15, weight='bold')
ax = plt.gca()
ax.patch.set_visible(False)
fig.patch.set_visible(False)
plt.show()
########################################################################################################################
########################################################################################################################
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
plt.style.use('classic')
dt_labels_tmstp = delta_g_hat.real[0]
rcParams['font.weight'] = 'bold'
dt = dt_labels_tmstp
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
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
plt.title("", fontsize=15, weight='bold')
plt.title("delta_g_hat", fontsize=15, weight='bold')
plt.xlabel("East [km]", fontsize=15, weight='bold')
plt.ylabel("North [km]", fontsize=15, weight='bold')
ax = plt.gca()
ax.patch.set_visible(False)
fig.patch.set_visible(False)
plt.show()
########################################################################################################################
########################################################################################################################
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
plt.style.use('classic')
t = delta_g_fft.real
dt_labels_tmstp = delta_g_2.real
rcParams['font.weight'] = 'bold'
dt = dt_labels_tmstp
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111)
ax.set_aspect('auto')
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
plt.title("", fontsize=15, weight='bold')
plt.title("delta_g", fontsize=15, weight='bold')
plt.xlabel("East [km]", fontsize=15, weight='bold')
plt.ylabel("North [km]", fontsize=15, weight='bold')
ax = plt.gca()
ax.patch.set_visible(False)
fig.patch.set_visible(False)
plt.show()
########################################################################################################################
########################################################################################################################
# end_time = time.time()
# execution_time = start_time - end_time
# print("Execution time:", execution_time)
