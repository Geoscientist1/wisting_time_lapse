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

app = dash.Dash(__name__)


########################################################################################################################

# This script models a simple model to be tested in SBLwiz. That is done in order to find out what the cause of
# Petrel models not being able to be read in SBLwiz.

########################################################################################################################

def simplemodelmaker():
    model = np.zeros((50, 50, 50))
    x = np.zeros((50, 50, 50))
    y = np.zeros((50, 50, 50))
    X = np.linspace(0, 5000, 50)
    Y = np.linspace(0, 5000, 50)
    z = np.zeros((50, 50, 50))

    dpt_increment = 0

    for k in range(0, 11):
        for j in range(50):
            for i in range(50):
                model[k, j, i] = 0.7
                z[k, j, i] = dpt_increment
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]
                if k == 10:
                    temp = z[k, j, i]
        dpt_increment -= 100

    dpt_increment = -100
    for k in range(11, 21):
        for j in range(50):
            for i in range(50):
                model[k, j, i] = 10
                z[k, j, i] = temp + dpt_increment
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]
                if k == 20:
                    temp1 = z[k, j, i]
        dpt_increment -= 100

    dpt_increment = -100
    for k in range(21, 37):
        for j in range(50):
            for i in range(50):
                model[k, j, i] = random.uniform(100, 100)
                z[k, j, i] = temp1 + dpt_increment
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]
                if k == 36:
                    temp2 = z[k, j, i]
        dpt_increment -= 100

    dpt_increment = -100
    for k in range(37, 50):
        for j in range(50):
            for i in range(50):
                model[k, j, i] = 20
                z[k, j, i] = temp2 + dpt_increment
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]
        dpt_increment -= 100

    return model, x, y, z


simple_model, x_simple, y_simple, z_simple = simplemodelmaker()


outF = open("myOutFile_" + "simple_model" + ".txt", "w")
for k in range(50):
    for j in range(50):
        for i in range(50):
            outF.write(str(i + 1) + " " + str(j + 1) + " " + str(k + 1) + " " + (str(x_simple[k, j, i]) +
                                                                                  " " + str(
                        y_simple[k, j, i]) + " " + str(
                        z_simple[k, j, i]) + " " + str(simple_model[k, j, i])))
            outF.write("\n")
outF.close()
