import random
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
import scipy.interpolate as si
import matplotlib.tri as mtri
import pandas as pd
import numpy as np
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import copy
from matplotlib import colors
from matplotlib import ticker
from matplotlib.patches import Rectangle
from numpy import *
import matplotlib.cm as cm
from matplotlib import cm
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
import plotly.graph_objects as go
import dash
# from atr_gen import atr_2727, atr_2732, atr_2737, atr_2742, atr_2747, atr_2752, atr_2757
# from dsw_dpt_slices import GRID, POR, dz, k_range, fluid_type, typo
from files_readme import phi, dpt, celltopdpt, f3, cel_cen_dpt

# from files_readme import A, B
A = 157
B = 159
C = 5

app = dash.Dash(__name__)

sat57 = np.loadtxt('sat57_coarse.txt', delimiter=' ')


########################################################################################################################

# This script reads a .txt file and substractes 1 from each and every index to have indices staring from 0 rather than 1
# This is done to see if that works with slbwiz forward modeling


# The function called distributor has as an objective to distribute data read from petrel
########################################################################################################################


def filemodifier(file):
    i = np.zeros(150401)  # The value of 150401 is pertaining to another model, the model of Claudia!
    j = np.zeros(150401)
    k = np.zeros(150401)
    x = np.zeros(150401)
    y = np.zeros(150401)
    z = np.zeros(150401)
    r = np.zeros(150401)
    line_count = 0
    index = 0
    with open(file, "r") as file:
        lines = file.readlines()
        file.close()
    for line in lines:
        arr = line.split()
        if line != "\n":
            line_count += 1
        i[index] = int(arr[0]) - 1
        j[index] = int(arr[1]) - 1
        k[index] = int(arr[2]) - 1
        x[index] = float(arr[3])
        y[index] = float(arr[4])
        z[index] = float(arr[5])
        r[index] = float(arr[6])
        index += 1
    return i, j, k, x, y, z, r


# file = 'hajar.txt'
# i, j, k, x, y, z, r = filemodifier(file)
#
# outF = open("myOutFile_" + "hajar_modified" + ".txt", "w")
# outF.write("\n")
# for index in range(150401):
#     outF.write(str(int(i[index])) + " " + str(int(j[index])) + " "
#                + str(int(k[index])) + " " + str(x[index]) + " " + str(y[index]) + " " + str(z[index]) + " " + str(
#         r[index]))
#     outF.write("\n")
# outF.close()


def top_reservoir_map():  # This script calculates and outputs the top reservoir map of the model
    # ( depth to top reservoir at that location in meters)
    M = A * B * 5  # Recall that M might need to be changed as a result of k=101 is actually k=k_range. That
    # can be fixed later.
    r = np.zeros((B, A))
    cel_dpt = np.reshape(cel_cen_dpt, (101, B, A))
    m = 0
    for k in range(1):
        for j in range(B):
            for i in range(A):
                if cel_dpt[k, j, i] != 0:
                    r[j, i] = cel_dpt[k, j, i]
                    m += 1
                else:
                    m += 1
    # rt = np.reshape(r, (B, A))
    # stk_top_res = np.zeros((B, A))
    # avg = round((M - x) / (A * B))
    # for i in range(A):
    #     for j in range(B):
    #         for k in range(101):
    #             stk_top_res[j, i] += rt[k, j, i]
    #         stk_top_res[j, i] /= avg
    result = r
    # for j in range(B):
    #     for i in range(A):
    #         if result[j, i] == 0:
    #             result[j, i] = 0.1
    return result


r = top_reservoir_map()


def thickness_reservoir_map(rt):  # This script calculates and outputs a map of the thickness of the
    # reservoir of the model (Thickness of reservoir above contact at each location in meters)
    dpt = np.reshape(cel_cen_dpt, (101, B, A))
    tck1 = np.zeros((B, A))
    tck2 = np.zeros((B, A))
    tck1_min = 717  # Mean Value
    tck2_max = 717  # Mean Value
    for j in range(B):
        for i in range(A):
            for k in range(1):
                if dpt[k, j, i] != 0:
                    variable_1 = dpt[k, j, i]
                    if variable_1 < tck1_min:
                        tck1[j, i] = variable_1
                        tck1_min = variable_1
                    else:
                        tck1[j, i] = tck1_min
                else:
                    tck1[j, i] = tck1_min
            tck1_min = 717  # Mean Value

    for j in range(B):
        for i in range(A):
            for k in range(8):
                if dpt[k, j, i] != 0:
                    variable_2 = dpt[k, j, i]
                    if variable_2 > tck2_max:
                        tck2[j, i] = variable_2
                        tck2_max = variable_2
                    else:
                        tck2[j, i] = tck2_max
                else:
                    tck2[j, i] = tck2_max
            tck2_max = 717  # Mean Value

    tck = tck2 - tck1
    for j in range(B):
        for i in range(A):
            if tck[j, i] > 37:
                tck[j, i] = 37
    # for k in range(1):
    #     for j in range(B):
    #         for i in range(A):
    #             W = celltopdpt[m]
    #             V = rt[j, i]
    #             # if celltopdpt[m] != 0:
    #             #     t[m] = celltopdpt[m] - rt[j, i]
    #             #     m += 1
    #             # else:
    #             #     t[m] = 0.1
    #             #     m += 1
    #             w = tck1[m]
    #             v = tck2[m]
    #             t[m] = tck2[m] - tck1[m]
    #             m += 1
    # th = np.reshape(t, (B, A))
    # for j in range(B):
    #     for i in range(A):
    #         if result[j, i] <= 0:
    #             result[j, i] = 0.1
    return tck


# tck = thickness_reservoir_map(r)

# X = np.linspace(598353.72, 611649.17, A)
# X = np.linspace(0, 11367.59, A)
# Y = np.linspace(8146346.56, 8160662.48, B)
# Y = np.linspace(0, 11557.04, B)


# for i in range(A):
#     X[i] =

# angle = -18.3
# results = [atr_2727, atr_2732, atr_2737, atr_2742, atr_2747, atr_2752, atr_2757]
# labels = ["atr_2727", "atr_2732", "atr_2737", "atr_2742", "atr_2747", "atr_2752", "atr_2757"]
# for i in range(len(results)):
#     dt = results[i]
#     outF0 = open(labels[i] + "_map" + ".txt", "w")
#     for j in range(B):
#         for i in range(A):
#             if dt[j, i] and r[j, i] == 0.1:
#                 S = r[j, i]
#                 continue
#             else:
#                 outF0.write("{0} {1} {2} {3} {4}".format(
#                     str(598353.72 - 1800 + (np.cos(np.deg2rad(angle)) * X[i] - np.sin(np.deg2rad(angle)) * Y[B - 1 - j])),
#                     str(8146346.56 + 3650 + (np.sin(np.deg2rad(angle)) * X[i] + np.cos(np.deg2rad(angle)) * Y[B - 1 - j])),
#                     str(round(dt[j, i])), str(round(r[j, i])), str(round(tck[j, i]))))
#
#             outF0.write("\n")
#     outF0.close()

# outF1 = open("myOutFile_" + "thickness_reservoir_wisting" + ".xyz", "w")
# for j in range(B):
#     for i in range(A):
#         if th[j, i] == 0.1:
#             continue
#         else:
#             outF1.write("{0} {1} {2}".format(str(X[i]), str(Y[j]), str(th[j, i])))
#
#         outF1.write("\n")
# outF1.close()


angle = -18.3
X = np.linspace(0, 11367.59, A)
Y = np.linspace(0, 11557.04, B)
Z = np.linspace(-596.02, -1082.86, C)
from Gravity_mod_grid_coaresening import sat57_coarse

dt = sat57_coarse
outF0 = open("sat57_coarse_coordinates" + ".txt", "w")
m = 0
extra = 0
for k in range(C):
    for i in range(A):
        for j in range(B):
            if dt[m] == 0.0:
                dt[m] = 1.0
                outF0.write("{0} {1} {2} {3}".format(
                    str(X[i]), str(Y[B - 1 - j]), str(r[j, i] + extra), str((dt[m]))))
                outF0.write("\n")
                m += 1
            else:
                outF0.write("{0} {1} {2} {3}".format(
                    str(X[i]), str(Y[B - 1 - j]), str(r[j, i] + extra), str((dt[m]))))
                outF0.write("\n")
                m += 1
    extra -= 10

outF0.close()
