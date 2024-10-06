import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import interpolate
from numpy import *
import cv2
import matplotlib.cm as cm
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from dsw_dpt_slices import k_range
from files_readme import A, B, cel_cen_dpt, dx_accum, dy_accum, phi, dx, dy, dz

rho_oil = 840  # Oil density
rho_gas = 150  # Gas density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
K = [388, 396, 404, 412, 420, 428, 436, 444, 452]  # 9 different vertical layers for the bathymetry ranging between


# 388 and 452 m with an increment equal to 8 meters

# The function calculates the Frechet kernel (Green's functions) for each single cell to be used in gravity inversion
def frechet_krn():
    XSCALE = 84
    YSCALE = 89
    length = len(K)
    n = 0
    a = G * (cel_cen_dpt - K[n]) * phi
    b = dx * dy * dz
    c = a * b
    frechet_krnl = np.zeros((B, A))
    for y_step in range(0, B):
        for x_step in range(0, A):
            print("y", y_step)
            print("x", x_step)
            m = 0
            for k in range(k_range):
                for j in range(B):
                    for i in range(A):
                        # print("h", h_step, "k", k, "j", j, "i", i)
                        r = np.sqrt(((x_step * XSCALE) - dx_accum[m]) ** 2 + ((y_step * YSCALE) - dy_accum[m]) ** 2)
                        d = pow((pow(r, 2) + pow(cel_cen_dpt[m] - K[n], 2)), 1.5)
                        frechet_krnl[158 - y_step, x_step] = ((c[m]) / d)
                        m += 1
    return frechet_krnl


frechet_kernel = frechet_krn()
for i in range(len(K)):
    np.savetxt('frechet_kernel_dpt_level_' + str(388) + '.txt', frechet_kernel, delimiter=' ',
               header='Frechet Kernel for the first bathymetry level: 388')
