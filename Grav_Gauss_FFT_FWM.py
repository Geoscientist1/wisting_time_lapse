import math
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import time
from scipy.ndimage import gaussian_filter
from scipy.ndimage import zoom
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from matplotlib import rcParams

start_time = time.time()
########################################################################################################################
rho_oil = 835  # Oil density
rho_water = 1040  # Water density
G = 6.67 * 1e-11  # Gravity constant
# Define the specific depth for observation surface
z0 = 0
########################################################################################################################
extra_term = 0
########################################################################################################################
M = 157  # This is the # of grid cells in the x-direction
x_min = 0
x_max = 13213.48
delta_x = (x_max - x_min) / M
xi = np.arange(0, (M + 1) * delta_x, delta_x)
########################################################################################################################
N = 159 # This is the # of grid cells in the y-direction
y_min = 0
y_max = 14173.40
delta_y = (y_max - y_min) / N
yi = np.arange(0, (N + 1) * delta_y, delta_y)
########################################################################################################################
L = 5  # This is the # of grid cells in the z-direction
z_min = -1082.86
z_max = -596.02
delta_z = (z_max - z_min) / L
zi = np.arange(0, (L + 1) * delta_z, delta_z)
########################################################################################################################
x_value = np.fft.fftfreq(M)
y_value = np.fft.fftfreq(N)
z_value = np.fft.fftfreq(L)
kx_array = zeros((L, N, M), dtype=float)
ky_array = zeros((L, N, M), dtype=float)
kz_array = zeros((L, N, M), dtype=float)
x_length = x_max - x_min
y_length = y_max - y_min
z_length = z_max - z_min
# now the loops to calculate the wavenumbers
for dim in range(L):
    for column in range(N):
        for row in range(M):
            kx_array[dim, column, row] = (2.0 * pi * x_value[row]) / (M * delta_x)
            ky_array[dim, column, row] = (2.0 * pi * y_value[column]) / (N * delta_y)
            kz_array[dim, column, row] = (2.0 * pi * z_value[dim]) / (L * delta_z)

# now for any row,column pair kx_array , and ky_array will hold the wavedomain coordinates
# of the correspoing point in some_data_wavedomain
########################################################################################################################
# Gaussian quadrature nodes and weights
n = 4
x, wi = np.polynomial.legendre.leggauss(n)
nodes = np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3, -1).T
wei = (wi * wi[:, None]).ravel()
weights = (wei * wi[:, None]).ravel()

########################################################################################################################
# Define the modified term function to handle division by zero
def modified_term(k, delta):
    epsilon = 1e-10
    return np.where(np.abs(k) < epsilon, delta / 2, np.sin(0.5 * k * delta) / k)


########################################################################################################################
def gauss_fft_fwm(delta_dsw_input, f):
    print("Gauss FFT Forward Modeling computation...")
    print("Inversion iteration: " + str(f) + " is running...")

    delta_rho_input = delta_dsw_input * (rho_water - rho_oil)
    delta_rho = np.reshape(delta_rho_input, (L, N, M))
    # Cutoff wavenumbers for kx and ky
    cutoff_wavenumber_kx = np.pi / delta_x  # Example value for kx cutoff
    cutoff_wavenumber_ky = np.pi / delta_y  # Example value for ky cutoff

    # Perform 2D FFT on each horizontal slice to get frequency domain representation
    rho_k_xy = np.fft.fft2(delta_rho, axes=(1, 2))

    # Generate wavenumber grids for the model (for x and y directions)
    kx = np.fft.fftfreq(M, d=delta_x) * 2 * np.pi
    ky = np.fft.fftfreq(N, d=delta_y) * 2 * np.pi

    # Meshgrid for wavenumbers in x and y directions
    KX, KY = np.meshgrid(ky, kx, indexing='ij')
    k_squared_x = KX ** 2
    k_squared_y = KY ** 2
    k_squared_xy = k_squared_x + k_squared_y

    # Apply different cutoff wavenumbers for kx and ky
    mask_x = k_squared_x <= cutoff_wavenumber_kx ** 2
    mask_y = k_squared_y <= cutoff_wavenumber_ky ** 2
    mask = mask_x & mask_y

    # Avoid very small values in k_squared to prevent division by zero
    k_squared_xy = np.clip(k_squared_xy, 1e-3, None)
    k_xy = np.sqrt(k_squared_xy)

    # Calculate modified terms for the model (for x and y directions)
    term_x = modified_term(KX, delta_x)
    term_y = modified_term(KY, delta_y)

    # Compute the Green's function term with correct sign convention for the model
    G_k_xy = (2 * np.pi * G / k_xy) * 2.0 * term_x * 2.0 * term_y

    # Initialize the array for storing results after applying the Green's function
    delta_g_k_xy = np.zeros_like(rho_k_xy[0, :, :], dtype=np.complex128)

    # Integrate the contributions from all z layers
    for k in range(L):
        # Apply Green's function for each z slice
        exp_sinh = np.exp(k_xy * (z0 - zi[k])) * 2.0 * np.sinh(k_xy * delta_z / 2.0) / delta_x / delta_y
        delta_g_k_xy_slice = mask * G_k_xy * rho_k_xy[k, :, :] * exp_sinh
        delta_g_k_xy += delta_g_k_xy_slice  # Sum the contributions

        # Define a Gaussian low-pass filter in the wavenumber domain
        alpha = 0.1  # Adjust alpha to control the filter sharpness
        H = np.exp(-alpha * k_squared_xy)

        # Apply the Gaussian low-pass filter
        delta_g_k_xy_filtered = delta_g_k_xy * H

        # Perform 2D IFFT to get the gravitational field anomaly in the spatial domain for each horizontal slice
        delta_g_filtered = np.fft.ifft2(delta_g_k_xy_filtered, axes=(0, 1)).real

        # Convert to µGal
        delta_g_z10_filtered_muGal = delta_g_filtered * 1e8  # Convert from m/s^2 to µGal

        # Optionally apply Gaussian filter to smooth the result
        Z_filtered = gaussian_filter(delta_g_z10_filtered_muGal, sigma=3, mode='reflect')

        f_result = Z_filtered / (4 * np.pi ** 2)

        # Making the output array and hence the measurement grid coarser
        zoom_factors = (57 / 159, 55 / 157)
        output_array = zoom(f_result, zoom_factors)

    return output_array


########################################################################################################################

# result = gauss_fft_fwm(input)
########################################################################################################################
# slice_index = 0
# X = np.linspace(0, 13213.48, 157)
# Y = np.linspace(0, 14173.40, 159)
# plt.style.use('classic')
# dt_labels_tmstp = result
# rcParams['font.weight'] = 'bold'
# dt = dt_labels_tmstp
# fig = plt.figure(figsize=(6.4, 4.8))
# ax = fig.add_subplot(111)
# ax.set_aspect('auto')
# Z = gaussian_filter(dt, sigma=0, mode='reflect')
# plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
#            cmap=plt.cm.get_cmap('jet', 1000),
#            interpolation='nearest', origin='upper')
# clb = plt.colorbar()
# plt.clim(0, 13)
# tick_locator = ticker.MaxNLocator(nbins=7)
# clb.locator = tick_locator
# clb.update_ticks()
# plt.title("", fontsize=15, weight='bold')
# plt.title("FFT Modeled data", fontsize=15, weight='bold')
# plt.xlabel("East [km]", fontsize=15, weight='bold')
# plt.ylabel("North [km]", fontsize=15, weight='bold')
# ax = plt.gca()
# ax.patch.set_visible(False)
# fig.patch.set_visible(False)
# plt.show()
# ########################################################################################################################
# X = np.linspace(0, 13213.48, 157)
# Y = np.linspace(0, 14173.40, 159)
# plt.style.use('classic')
# dt_labels_tmstp = delta_g
# rcParams['font.weight'] = 'bold'
# dt = dt_labels_tmstp
# fig = plt.figure(figsize=(6.4, 4.8))
# ax = fig.add_subplot(111)
# ax.set_aspect('auto')
# Z = gaussian_filter(dt, sigma=0, mode='reflect')
# plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
#            cmap=plt.cm.get_cmap('jet', 1000),
#            interpolation='nearest', origin='upper')
# clb = plt.colorbar()
# plt.clim(0, 13)
# tick_locator = ticker.MaxNLocator(nbins=7)
# clb.locator = tick_locator
# clb.update_ticks()
# clb.ax.set_title('$(\u03BC$Gal)', fontsize=15, weight='bold')
# plt.title("", fontsize=15, weight='bold')
# plt.title("Observed data", fontsize=15, weight='bold')
# plt.xlabel("East [km]", fontsize=15, weight='bold')
# plt.ylabel("North [km]", fontsize=15, weight='bold')
# ax = plt.gca()
# ax.patch.set_visible(False)
# fig.patch.set_visible(False)
# plt.show()
# ########################################################################################################################
# X = np.linspace(0, 13213.48, 157)
# Y = np.linspace(0, 14173.40, 159)
# plt.style.use('classic')
# dt_labels_tmstp = np.absolute(result - delta_g)
# rcParams['font.weight'] = 'bold'
# dt = dt_labels_tmstp
# fig = plt.figure(figsize=(6.4, 4.8))
# ax = fig.add_subplot(111)
# ax.set_aspect('auto')
# Z = gaussian_filter(dt, sigma=3, mode='reflect')
# plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
#            cmap=plt.cm.get_cmap('jet', 1000),
#            interpolation='nearest', origin='upper')
# clb = plt.colorbar()
# plt.clim(0, 4)
# tick_locator = ticker.MaxNLocator(nbins=7)
# clb.locator = tick_locator
# clb.update_ticks()
# clb.ax.set_title('$(\u03BC$Gal)', fontsize=15, weight='bold')
# plt.title("", fontsize=15, weight='bold')
# plt.title("Difference", fontsize=15, weight='bold')
# plt.xlabel("East [km]", fontsize=15, weight='bold')
# plt.ylabel("North [km]", fontsize=15, weight='bold')
# ax = plt.gca()
# ax.patch.set_visible(False)
# fig.patch.set_visible(False)
# plt.show()
# ########################################################################################################################
# #######################################################################################################################
# # end_time = time.time()
# # execution_time = start_time - end_time
# # print("Execution time:", execution_time)
