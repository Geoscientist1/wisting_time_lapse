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

# This script reads an .earthVision grid file from Petrel and get rids of the row and column indices, such that we have
# a pure xyz file to work with.
# That will allow the file to be used in Brickville for model building.


# The function called deleter has as an objective to delete row and column indices
########################################################################################################################
def deleter(file):
    string_filename = file
    with open(file, "r") as file:
        lines = file.readlines()
        file.close()
    outF = open("wow.xyz", "w")
    for line in lines:
        arr = line.split()
        outF.write(str(arr[0]) + " " + str(arr[1]) + " " + str(arr[2]))
        outF.write("\n")
    outF.close()


file = 'wow.xyz'
deleter(file)
