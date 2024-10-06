# import numpy as np
# N = 7
# X = np.linspace(595353.72, 611949.17, N)
# Y = np.linspace(8145346.56, 8163662.48, N)
#
# outF0 = open("MCPL_RxPos" + ".xyz", "w")
# for i in range(N):
#     for j in range(N):
#         outF0.write("{0} {1} {2}".format(str(round(X[i])), str(round(Y[j])), str(0.0)))
#         outF0.write("\n")
# outF0.close()
#
# frechet_kernel = np.loadtxt("frechet_kernel_inc_10_dpt_slice7.txt", delimiter=" ")
# print("This is the length of the file:", len(frechet_kernel))
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import cv2
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from dsw_dpt_slices import GRID, POR, dz, k_range, fluid_type, typo
from files_readme import celltopdpt, cell_vol
from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
from files_readme import dx, dy, dz, A, B, C, phi, porv, cel_cen_dpt, dx, dy, dz, dx_accum, dy_accum

rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
dv = dx * dy * dz
Phi = phi
Porv = porv

# Now we account for the different water saturations corresponding to the different times from year 2027 all the way
# to year 2039
################################################################################
################################################################################
frechet_kernel = np.loadtxt("frechet_kernel_dpt_level_388.txt", delimiter=" ")


def gravity_mod(sat_1, sat_2, soil_1, soil_2, frechet_kernel):
    rho_1 = ((1 - soil_1) * rho_water) + ((1 - sat_1) * rho_oil)
    rho_2 = ((1 - soil_2) * rho_water) + ((1 - sat_2) * rho_oil)
    delta_rho = rho_2 - rho_1
    delta_rho = np.reshape(delta_rho, (C, B, A))
    # frechet_kernel = np.reshape(frechet_kernel, (k_range, B, A))
    delta_g = np.zeros((B, A))
    dV = np.reshape(dv, (C, B, A))
    phi = np.reshape(Phi, (C, B, A))
    for i in range(A):
        for j in range(B):
            for k in range(k_range):
                print("j", j)
                print("i", i)
                # P = frechet_kernel[k, j, i]
                Q = delta_rho[k, j, i]
                R = dV[k, j, i]
                T = phi[k, j, i]
                # S = P * Q * R * T
                delta_g[j, i] += ((frechet_kernel[j, i]) * delta_rho[k, j, i] * dV[k, j, i] * phi[k, j, i])
    return delta_g * 1e8


################################################################################
print('Nooooooooooooooooooooooooooooooooooooooooooooooooooow')
gv27_27 = gravity_mod(sat27, sat27, soil_27, soil_27, frechet_kernel)
gv27_29 = gravity_mod(sat27, sat29, soil_27, soil_29, frechet_kernel)
gv27_31 = gravity_mod(sat27, sat31, soil_27, soil_31, frechet_kernel)
gv27_33 = gravity_mod(sat27, sat33, soil_27, soil_33, frechet_kernel)
gv27_35 = gravity_mod(sat27, sat35, soil_27, soil_35, frechet_kernel)
gv27_39 = gravity_mod(sat27, sat39, soil_27, soil_39, frechet_kernel)
gv27_45 = gravity_mod(sat27, sat45, soil_27, soil_45, frechet_kernel)
gv27_49 = gravity_mod(sat27, sat49, soil_27, soil_49, frechet_kernel)
gv27_53 = gravity_mod(sat27, sat53, soil_27, soil_53, frechet_kernel)
gv27_57 = gravity_mod(sat27, sat57, soil_27, soil_57, frechet_kernel)


np.savetxt('grav_data_27_frechet.txt', gv27_27, delimiter=' ', header='Gravity data for year 2027')
np.savetxt('grav_data_29_frechet.txt', gv27_29, delimiter=' ', header='Gravity data for year 2029')
np.savetxt('grav_data_31_frechet.txt', gv27_31, delimiter=' ', header='Gravity data for year 2031')
np.savetxt('grav_data_33_frechet.txt', gv27_33, delimiter=' ', header='Gravity data for year 2033')
np.savetxt('grav_data_35_frechet.txt', gv27_35, delimiter=' ', header='Gravity data for year 2035')
np.savetxt('grav_data_39_frechet.txt', gv27_39, delimiter=' ', header='Gravity data for year 2039')
np.savetxt('grav_data_45_frechet.txt', gv27_45, delimiter=' ', header='Gravity data for year 2045')
np.savetxt('grav_data_49_frechet.txt', gv27_49, delimiter=' ', header='Gravity data for year 2049')
np.savetxt('grav_data_53_frechet.txt', gv27_53, delimiter=' ', header='Gravity data for year 2053')
np.savetxt('grav_data_57_frechet.txt', gv27_57, delimiter=' ', header='Gravity data for year 2057')

gv_f_27 = np.loadtxt("grav_data_27_frechet.txt", delimiter=" ")
gv_f_29 = np.loadtxt("grav_data_29_frechet.txt", delimiter=" ")
gv_f_31 = np.loadtxt("grav_data_31_frechet.txt", delimiter=" ")
gv_f_33 = np.loadtxt("grav_data_33_frechet.txt", delimiter=" ")
gv_f_35 = np.loadtxt("grav_data_35_frechet.txt", delimiter=" ")
gv_f_39 = np.loadtxt("grav_data_39_frechet.txt", delimiter=" ")
gv_f_45 = np.loadtxt("grav_data_45_frechet.txt", delimiter=" ")
gv_f_49 = np.loadtxt("grav_data_49_frechet.txt", delimiter=" ")
gv_f_53 = np.loadtxt("grav_data_53_frechet.txt", delimiter=" ")
gv_f_57 = np.loadtxt("grav_data_57_frechet.txt", delimiter=" ")
plt.style.use('classic')
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
XSCALE = 84
YSCALE = 89

dt_labels_tmstp = [gv_f_27, gv_f_29, gv_f_31, gv_f_33, gv_f_35, gv_f_39, gv_f_45, gv_f_49, gv_f_53, gv_f_57]
for i in range(10):
    dt = dt_labels_tmstp[i]
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')
    f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
    Z = gaussian_filter(dt, sigma=0, mode='reflect')
    plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
               cmap=plt.cm.get_cmap('jet', 1000),
               interpolation='nearest', origin='upper')
    clb = plt.colorbar()
    tick_locator = ticker.MaxNLocator(nbins=7)
    clb.locator = tick_locator
    clb.update_ticks()
    clb.ax.set_title('$(\u03BC$Gal)', fontsize=15)
    ax.set_ylabel("Y(m)", labelpad=15)
    plt.clim(0, 10)
    plt.xlabel("X [$km$]")
    plt.ylabel("Y [$km$]")
    ax = plt.gca()
    ax.patch.set_visible(False)
    fig.patch.set_visible(False)
    ax.grid(color='b', linestyle='-', linewidth=0.5)
    plt.show()
