import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import interpolate
from numpy import *
# import cv2
import matplotlib.cm as cm
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from files_readme import dx, dy, dz, dy_avg, dz_avg, phi
from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
from files_readme import sgas_27, sgas_29, sgas_31, sgas_33, sgas_35, sgas_39, sgas_45, sgas_49, sgas_53, sgas_57
from files_readme import A, B

rho_oil = 840  # Oil density
rho_gas = 0.75  # Gas density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
k_range = 9  # Number of cells considered along the z-direction, GOC is set at cell 21, and OWC is set at cell 51


def check(list1, val):
    for x in list1:
        if val == x:
            return False
    return True


# The function calculating water saturation from density changes
def grav_rev(sat_1, sat_2):
    N = len(sat_1)
    rho_1 = ((1 - sat_1) * rho_oil) + (sat_1 * rho_water)
    rho_2 = ((1 - sat_2) * rho_oil) + (sat_2 * rho_water)
    delta_rho = rho_2 - rho_1
    dens = np.zeros((B, A))
    S = np.zeros(N)
    m = 0

    for k in range(101):
        for j in range(B):
            for i in range(A):
                if phi[m] != 0:
                    temp = delta_rho[m] / ((rho_water - rho_oil) * phi[m])
                    dens[j, i] = 1 - (sum(temp) / 101)
                    S[m] = temp
                    m += 1
                else:
                    S[m] = -1
                    m += 1

    Sw = np.reshape(S, (101, B, A))
    dSw = np.zeros((B, A))

    for i in range(A):
        for j in range(B):
            for n in range(101):
                dSw[j, i] += Sw[n, j, i]
            dSw[j, i] = dSw[j, i] / 101

    return dens, dSw


########################################################################################################################
def grid(d_z, phi):
    if k_range != 101:
        m = 0
        dz = np.zeros(k_range * A * B)
        for k in range(k_range):
            for j in range(B):
                for i in range(A):
                    dz[m] = d_z[m]
                    m += 1
        dz = np.reshape(dz, (k_range, B, A))
    else:
        dz = np.reshape(d_z, (k_range, B, A))

    dZ = np.zeros((10, B, A))
    phi = np.reshape(phi, (101, B, A))
    poro = np.zeros((10, B, A))

    # s = 0
    # fk = 1
    # temp0 = 0
    # temp1 = 0
    # for i in range(157):
    #     for j in range(159):
    #         for n in range(k_range):
    #             temp0 += dz[n, j, i]
    #             temp1 += phi[n, j, i]
    #             if n != 0 and mod(n, 9) == 0 and n != 99:
    #                 dZ[s, j, i] = temp0 / (fk * 10)
    #                 poro[s, j, i] = temp1 / (fk * 10)
    #                 s += 1
    #                 fk += 1
    #                 temp0 = 0
    #                 temp1 = 0
    #         s = 0
    #         fk = 1
    return dz, poro


########################################################################################################################

def dsw_slices(fluid_saturation, sat_1, sat_2, soil_1, soil_2, sgas_1, sgas_2):
    N = len(sat_1)
    delta_sw = sat_2 - sat_1
    rho_avg_w = rho_water - (rho_oil + rho_gas)
    delta_so = soil_2 - soil_1
    rho_avg_o = rho_oil - (rho_water + rho_gas)
    delta_sg = sgas_2 - sgas_1
    rho_avg_g = rho_gas - (rho_water + rho_oil)

    if fluid_saturation == "Water-Oil_saturation":
        rho_1 = ((1 - soil_1) * rho_water) + ((1 - sat_1) * rho_oil)
        rho_2 = ((1 - soil_2) * rho_water) + ((1 - sat_2) * rho_oil)
        delta_rho = rho_2 - rho_1

    elif fluid_saturation == "Water-Gas_saturation":
        rho_1 = ((1 - sgas_1) * rho_water) + ((1 - sat_1) * rho_gas)
        rho_2 = ((1 - sgas_2) * rho_water) + ((1 - sat_2) * rho_gas)
        delta_rho = rho_2 - rho_1

    #     rho_1 = (sat_1 * rho_water) + ((1 - sat_1) * rho_oil)
    #     rho_2 = (sat_2 * rho_water) + ((1 - sat_2) * rho_oil)
    #     delta_rho = rho_2 - rho_1
    # elif fluid_saturation == "Oil_saturation":
    #     rho_1 = (soil_1 * rho_oil) + ((1 - soil_1) * rho_water)
    #     rho_2 = (soil_2 * rho_oil) + ((1 - soil_2) * rho_water)
    #     delta_rho = rho_1 - rho_2
    # elif fluid_saturation == "Gas_saturation":

    m = 0
    S = np.zeros(k_range * A * B)

    for k in range(k_range):
        for j in range(B):
            for i in range(A):
                if phi[m] != 0:
                    if fluid_saturation == "Water_saturation":
                        S[m] = (sat_2[m] - sat_1[m])
                        m += 1
                    elif fluid_saturation == "Oil_saturation":
                        S[m] = np.abs(soil_2[m] - soil_1[m]) / phi[m]
                        m += 1
                    elif fluid_saturation == "Gas_saturation":
                        S[m] = np.abs(sgas_2[m] - sgas_1[m])
                        m += 1
                    elif fluid_saturation == "Water-Oil_saturation":
                        S[m] = delta_rho[m] / ((rho_water - rho_oil) * phi[m])
                        m += 1
                    elif fluid_saturation == "Water-Gas_saturation":
                        S[m] = delta_rho[m] / ((rho_water - rho_gas) * phi[m])
                        # if S[m] < 0:
                        #     S[m] = 0
                        m += 1
                else:
                    S[m] = -1
                    m += 1

    Sw = np.reshape(S, (k_range, B, A))
    dsw = np.zeros((B, A))
    dSw = np.zeros((10, B, A))
    s = 0
    fk = 1
    temp = 0
    # for i in range(157):
    #     for j in range(159):
    #         for n in range(k_range):
    #             temp += Sw[n, j, i]
    #             if n != 0 and mod(n, 5) == 0 and n != 99:
    #                 dSw[s, j, i] = temp / (fk * 5)
    #                 s += 1
    #                 temp = 0  # inside the loop or outside depends on whether we need the accumulated average or
    #         # single layer average
    #         s = 0
    #         fk = 1
    return Sw  # return gas, oil or water saturation depending on the fluid_type "fluid_saturation"


########################################################################################################################

fluid_type = "Water_saturation"
typo = "stk"

########################################################################################################################
GRID, POR = grid(dz, phi)

dsw_2727 = dsw_slices(fluid_type, sat27, sat27, soil_27, soil_27, sgas_27, sgas_27)
dsw_2729 = dsw_slices(fluid_type, sat27, sat29, soil_27, soil_29, sgas_27, sgas_29)
dsw_2731 = dsw_slices(fluid_type, sat27, sat31, soil_27, soil_31, sgas_27, sgas_31)
dsw_2733 = dsw_slices(fluid_type, sat27, sat33, soil_27, soil_33, sgas_27, sgas_33)
dsw_2735 = dsw_slices(fluid_type, sat27, sat35, soil_27, soil_35, sgas_27, sgas_35)
dsw_2739 = dsw_slices(fluid_type, sat27, sat39, soil_27, soil_39, sgas_27, sgas_39)
dsw_2745 = dsw_slices(fluid_type, sat27, sat45, soil_27, soil_45, sgas_27, sgas_45)
dsw_2749 = dsw_slices(fluid_type, sat27, sat49, soil_27, soil_49, sgas_27, sgas_49)
dsw_2753 = dsw_slices(fluid_type, sat27, sat53, soil_27, soil_53, sgas_27, sgas_53)
dsw_2757 = dsw_slices(fluid_type, sat27, sat57, soil_27, soil_57, sgas_27, sgas_57)

# X = np.zeros(A)
# Y = np.zeros(B)
# for i in range(A):
#     for j in range(B):
#         for n in range(101):
#             dX[i, j] += d_x[n, j, i]
#             u += d_x[n, j, i]
#             dY[i, j] += d_y[n, j, i]
#             v += d_y[n, j, i]
#         dX[i, j] = (dX[i, j] / 101)
#         dY[i, j] = (dY[i, j] / 101)


########################################################################################################################

plt.style.use('classic')

X = np.linspace(0, 13296, A)
Y = np.linspace(0, 14316, B)
########################################################################################################################

if fluid_type == "Gas_saturation":

    oupt_labels_tmstp = ["dsg_2727", "dsg_2729", "dsg_2731", "dsg_2733", "dsg_2735", "dsg_2739", "dsg_2745", "dsg_2749",
                         "dsg_2753", "dsg_2757"]

elif fluid_type == "Water_saturation":

    oupt_labels_tmstp = ["dsw_2727", "dsw_2729", "dsw_2731", "dsw_2733", "dsw_2735", "dsw_2739", "dsw_2745", "dsw_2749",
                         "dsw_2753", "dsw_2757"]

elif fluid_type == "Oil_saturation":

    oupt_labels_tmstp = ["dso_2727", "dso_2729", "dso_2731", "dso_2733", "dso_2735", "dso_2739", "dso_2745", "dso_2749",
                         "dso_2753", "dso_2757"]

elif fluid_type == "Water-Oil_saturation":

    oupt_labels_tmstp = ["dswo_2727", "dswo_2729", "dswo_2731", "dswo_2733", "dwso_2735", "dswo_2739", "dswo_2745",
                         "dswo_2749", "dswo_2753", "dswo_2757"]

elif fluid_type == "Water-Gas_saturation":

    oupt_labels_tmstp = ["dswg_2727", "dswg_2729", "dswg_2731", "dswg_2733", "dwsg_2735", "dswg_2739", "dswg_2745",
                         "dswg_2749", "dswg_2753", "dswg_2757"]

########################################################################################################################
fig_labels_tmstp = ["27-27", "27-29", "27-31", "27-33", "27-35", "27-39", "27-45", "27-49", "27-53",
                    "27-57"]
dt_labels_tmstp = [dsw_2727, dsw_2729, dsw_2731, dsw_2733, dsw_2735, dsw_2739, dsw_2745, dsw_2749, dsw_2753,
                   dsw_2757]
fig_labels_layer = ["Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "Layer 7", "Layer 8", "Layer 9",
                    "Layer 10"]
oupt_labels_layer = ["_1", "_2", "_3", "_4", "_5", "_6", "_7", "_8", "_9", "_10"]
########################################################################################################################
# for i in range(10):
#     fig = plt.figure(figsize=(6.4, 4.8))
#     ax = fig.add_subplot(111)
#     ax.set_aspect('equal')
#     gr = GRID[i]
#     f = RectBivariateSpline(Y, X, gr)
#     Z = f(Y, X)
#     Z = gaussian_filter(gr, sigma=3)
#     plt.imshow(gr, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 75), interpolation='nearest',
#                origin='upper')
#     plt.colorbar(extend='max')
#     plt.clim(0, 5)
#     plt.title('dz--->' + fig_labels_layer[i])
#     plt.xlabel("X(m)")
#     ax.xaxis.tick_top()
#     plt.ylabel("Y(m)")
#     plt.show()
#     fig.savefig("dz" + oupt_labels_layer[i])
########################################################################################################################
# for i in range(10):
#     fig = plt.figure(figsize=(6.4, 4.8))
#     ax = fig.add_subplot(111)
#     ax.set_aspect('equal')
#     porosity = POR[i]
#     f = RectBivariateSpline(Y, X, porosity)
#     Z = f(Y, X)
#     Z = gaussian_filter(porosity, sigma=3)
#     plt.imshow(porosity, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 75),
#                interpolation='nearest', origin='upper')
#     plt.colorbar(extend='max')
#     plt.clim(0, 0.5)
#     plt.title('phi--->' + fig_labels_layer[i])
#     plt.xlabel("X(m)")
#     ax.xaxis.tick_top()
#     plt.ylabel("Y(m)")
#     plt.show()
#     fig.savefig("phi" + oupt_labels_layer[i])
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

# for i in range(10):
#     dt = dt_labels_tmstp[i]
#     for j in range(10):
#         fig = plt.figure(figsize=(6.4, 4.8))
#         ax = fig.add_subplot(111)
#         ax.set_aspect('equal')
#         f = RectBivariateSpline(Y, X, dt[j] / 10)
#         Z = f(Y, X)
#         Z = gaussian_filter(dt[j] / 10, sigma=3)
#         k = [-0.15, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0]
#
#         if fluid_type == "Gas_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('winter', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dsg ---> ' + fig_labels_tmstp[i] + ' ---> ' + fig_labels_layer[j] + "_" + phs)
#             plt.clim(k[j], 0.3)
#
#         elif fluid_type == "Water_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dsw ---> ' + fig_labels_tmstp[i] + ' ---> ' + fig_labels_layer[j] + "_" + phs)
#             plt.clim(-0.1, 0.3)
#
#         elif fluid_type == "Oil_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('winter', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dso ---> ' + fig_labels_tmstp[i] + ' ---> ' + fig_labels_layer[j] + "_" + phs)
#             plt.clim(k[j], 0.5)
#
#         elif fluid_type == "Water-Oil_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dswo ---> ' + fig_labels_tmstp[i] + ' ---> ' + fig_labels_layer[j] + "_" + phs)
#             plt.clim(-0.1, 0.3)
#
#         elif fluid_type == "Water-Gas_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dswg ---> ' + fig_labels_tmstp[i] + ' ---> ' + fig_labels_layer[j] + "_" + phs)
#             plt.clim(-0.3, 0.3)
#
#         plt.colorbar(extend='max')
#
#         plt.xlabel("X(m)")
#         ax.xaxis.tick_top()
#         plt.ylabel("Y(m)")
#         plt.show()
#         if fluid_type == "Water-Oil_saturation" or fluid_type == "Water-Gas_saturation":
#             fig.savefig(oupt_labels_tmstp[i] + oupt_labels_layer[j] + "_" + phs)
#         else:
#             fig.savefig(oupt_labels_tmstp[i] + oupt_labels_layer[j])
########################################################################################################################


# for i in range(10):
#     dt = dt_labels_tmstp[i]
#     for j in range(k_range):
#         fig = plt.figure(figsize=(6.4, 4.8))
#         ax = fig.add_subplot(111)
#         ax.set_aspect('equal')
#         f = RectBivariateSpline(Y, X, dt[j])
#         Z = f(Y, X)
#         Z = gaussian_filter(dt[j], sigma=3)
#         k = [-0.15, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0]
#
#         if fluid_type == "Gas_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('winter', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dsg ---> ' + fig_labels_tmstp[i] + ' ---> ' + "_layer " + str(j + 1) + "_" + phs)
#             plt.clim(k[j], 0.3)
#
#         elif fluid_type == "Water_saturation":
#             plt.imshow(dt[j], extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dsw ---> ' + fig_labels_tmstp[i] + ' ---> ' + "_layer " + str(j + 1) + "_" + phs)
#             plt.clim(-0.1, 0.75)
#
#         elif fluid_type == "Oil_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('winter', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dso ---> ' + fig_labels_tmstp[i] + ' ---> ' + "_layer " + str(j + 1) + "_" + phs)
#             plt.clim(k[j], 0.5)
#
#         elif fluid_type == "Water-Oil_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dswo ---> ' + fig_labels_tmstp[i] + ' ---> ' + "_layer " + str(j + 1) + "_" + phs)
#             plt.clim(-0.1, 0.3)
#
#         elif fluid_type == "Water-Gas_saturation":
#             plt.imshow(dt[j] / 10, extent=[min(X), max(X), max(Y), min(Y)], cmap=plt.cm.get_cmap('jet', 1000),
#                        interpolation='nearest', origin='upper')
#             plt.title('dswg ---> ' + fig_labels_tmstp[i] + ' ---> ' + "_layer " + str(j + 1) + "_" + phs)
#             plt.clim(-0.3, 0.3)
#
#         plt.colorbar(extend='max')
#
#         plt.xlabel("X(m)")
#         ax.xaxis.tick_top()
#         plt.ylabel("Y(m)")
#         plt.show()
#         if fluid_type == "Water-Oil_saturation" or fluid_type == "Water-Gas_saturation":
#             fig.savefig(oupt_labels_tmstp[i] + "__" + str(j + 1) + "_" + phs)
#         else:
#             fig.savefig(oupt_labels_tmstp[i] + "__" + str(j + 1))
