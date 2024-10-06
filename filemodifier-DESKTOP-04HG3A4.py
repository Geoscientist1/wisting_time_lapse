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
from atr_gen import atr_2727
from dsw_dpt_slices import GRID, POR, dz, k_range, fluid_type, typo
from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
from files_readme import sgas_27, sgas_29, sgas_31, sgas_33, sgas_35, sgas_39, sgas_45, sgas_49, sgas_53, sgas_57
from files_readme import phi, dpt, celltopdpt, f3

app = dash.Dash(__name__)


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


file = 'hajar.txt'
i, j, k, x, y, z, r = filemodifier(file)

outF = open("myOutFile_" + "hajar_modified" + ".txt", "w")
outF.write("\n")
for index in range(150401):
    outF.write(str(int(i[index])) + " " + str(int(j[index])) + " "
               + str(int(k[index])) + " " + str(x[index]) + " " + str(y[index]) + " " + str(z[index]) + " " + str(
        r[index]))
    outF.write("\n")
outF.close()


def top_reservoir_map(sat_1, sat_2):  # This script calculates and outputs the top reservoir map of the model
    # ( depth to top reservoir at that location in meters)
    M = 157 * 159 * 1  # Recall that M might need to be changed as a result of k=101 is actually k=k_range. That
    # can be fixed later.
    r = np.zeros(M)
    m = 0
    x = 0
    for k in range(1):
        for j in range(159):
            for i in range(157):
                if celltopdpt[m] != 0:
                    r[m] = celltopdpt[m]
                    m += 1
                else:
                    r[m] = 0.1
                    m += 1
    rt = np.reshape(r, (159, 157))
    # stk_top_res = np.zeros((159, 157))
    # avg = round((M - x) / (157 * 159))
    # for i in range(157):
    #     for j in range(159):
    #         for k in range(101):
    #             stk_top_res[j, i] += rt[k, j, i]
    #         stk_top_res[j, i] /= avg
    result = rt
    # for j in range(159):
    #     for i in range(157):
    #         if result[j, i] == 0:
    #             result[j, i] = 0.1
    return result


r = top_reservoir_map(sat27, sat57)


def thickness_reservoir_map(sat_1, sat_2, rt):  # This script calculates and outputs a map of the thickness of the
    # reservoir of the model (Thickness of reservoir above contact at each location in meters)
    M = 157 * 159 * 101  # Recall that M might need to be changed as a result of k=101 is actually k=k_range. That
    # can be fixed later.
    t = np.zeros(157 * 159)
    m = 0
    x = 0
    for k in range(101):
        for j in range(157):
            for i in range(159):
                if sat_1[m] != 0 and sat_2[m] != 0 and phi[m] != 0:
                    if sat_2[m] >= 0.8:
                        t[x] = celltopdpt[m]
                        m += 1
                        x += 1
                        break
                else:
                    m += 1
    th = np.reshape(t, (159, 157))
    result = th
    for j in range(159):
        for i in range(157):
            if result[j, i] <= 0:
                result[j, i] = 0.1
    return result


th = thickness_reservoir_map(sat27, sat27, r)

X = np.linspace(598353.72, 611649.17, 157)
Y = np.linspace(8146346.56, 8160662.48, 159)

outF0 = open("myOutFile_" + "combined" + ".txt", "w")
for j in range(159):
    for i in range(157):
        if atr_2727[j, i] and r[j, i] == 0.1:
            S = r[j, i]
            continue
        else:
            outF0.write("{0} {1} {2} {3} {4}".format(str(round(X[i])), str(round(Y[158 - j])),
                                                     str(round(atr_2727[j, i])), str(round(r[j, i])),
                                                     str(abs(round(r[j, i] - 697)))))

        outF0.write("\n")
outF0.close()

outF1 = open("myOutFile_" + "thickness_reservoir_wisting" + ".xyz", "w")
for j in range(159):
    for i in range(157):
        if th[j, i] == 0.1:
            continue
        else:
            outF1.write("{0} {1} {2}".format(str(X[i]), str(Y[j]), str(th[j, i])))

        outF1.write("\n")
outF1.close()
