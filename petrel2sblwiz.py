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

# This script models the background resistivity and imprints the gslbib imported petrel resistivity models into this
# background model. That is done in order to provide the petrel simulation grid models with some background
# resistivities above and below it to be used in csem forward modeling in sblwiz.

########################################################################################################################


def petrel2sblwiz(file, new_file):
    line_count = 0
    imported_model = np.zeros((101, 159, 157))

    with open(file, "r") as file:
        lines = file.readlines()
        file.close()

    del lines[0]
    del lines[0]
    del lines[0]
    del lines[0]
    del lines[0]
    del lines[0]
    del lines[0]
    del lines[0]
    del lines[0]

    # for line in lines:
    #     arr = line.split()
    #     if float(arr[6]) == 100:
    #         del lines[0]

    new_file = open(new_file, "w+")
    for line in lines:
        new_file.write(line)
        arr = line.split()
        if line != "\n":
            line_count += 1
        index = 0
        for s in range(len(arr)):

            i = int(arr[0]) - 1
            index += 1
            j = int(arr[1]) - 1
            index += 1
            k = 101 - int(arr[2])
            index += 1
            imported_model[k, j, i] = float(arr[6])
            if index == 3:
                index = 0
                continue

    new_file.close()

    return imported_model


########################################################################################################################


def multiple3dsurfaceplots(file):
    line_count = 0
    M = 518
    X = np.linspace(0, 13295.46, 157)
    Y = np.linspace(0, 14315.93, 159)
    model = np.zeros((M, 159, 157))
    x = np.zeros((M, 159, 157))
    y = np.zeros((M, 159, 157))
    z = np.zeros((M, 159, 157))
    dpt_increment = 0
    temp = 0
    temp1 = 0
    temp2 = 0

    # 134 cells for the water layer with uniform resistivity being equal to 0.1 ohm-m
    for k in range(0, 134):
        for j in range(159):
            for i in range(157):
                model[k, j, i] = 0.1
                z[k, j, i] = dpt_increment
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]
                if k == 133:
                    temp = z[k, j, i]
        dpt_increment -= 3

    # 83 cells for the overburden layer with uniform resistivity being equal to 10 ohm-m
    dpt_increment = -3
    for k in range(134, 217):
        for j in range(159):
            for i in range(157):
                model[k, j, i] = 10
                z[k, j, i] = temp + dpt_increment
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]
                if k == 216:
                    temp1 = np.amin(z[k, :, :])
        dpt_increment -= 3

    # 101 cells for the reservoir, where values are exported from a .gslib file from Petrel
    dpt_increment = -3
    dpt_increment_counter = 101
    k_control = 0
    with open(file, "r") as file:
        lines = file.readlines()
        file.close()
    for line in lines:
        arr = line.split()
        if line != "\n":
            line_count += 1
        index = 0
        i = int(arr[0]) - 1
        index += 1
        j = int(arr[1]) - 1
        index += 1
        k = 101 - int(arr[2]) + 216
        index += 1
        if float(arr[3]) and float(arr[4]) and float(arr[5]) and float(arr[6]) == -99999:
            model[k, j, i] = 15
            x[k, j, i] = X[i]
            y[k, j, i] = Y[j]
            z[k, j, i] = temp1 + dpt_increment
            if i == 156 and j == 158:
                dpt_increment -= 20
        else:
            model[k, j, i] = float(arr[6])
            s = float(arr[6])
            t = float(arr[5])
            o = temp1 + float(arr[5])
            x[k, j, i] = X[i]
            y[k, j, i] = Y[j]
            z[k, j, i] = temp1 + dpt_increment
            if i == 156 and j == 158:
                dpt_increment -= 20
        if k == 149:
            temp2 = np.amin(z[k, :, :])
        if index == 3:
            index = 0
            continue

    # 50 cells for the halfspace layer with uniform resistivity being equal to 10 ohm-m
    dpt_increment = -20
    for k in range(150, M):
        for j in range(159):
            for i in range(157):
                model[k, j, i] = 20
                z[k, j, i] = temp2 + dpt_increment
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]
        dpt_increment -= 20

    file.close()

    return model, x, y, z


########################################################################################################################
resis_target = 'moon.txt'
resis_target_edited = 'moon_edited.txt'
converted_model = petrel2sblwiz(resis_target, resis_target_edited)
plotted_model, x_coordinates, y_coordinates, z_coordinates = multiple3dsurfaceplots(resis_target_edited)
# handles = [Rectangle((0, 0), 1, 1, color=Color, ec="k")]

plt.style.use('classic')

X = np.linspace(598982.08, 610735.23, 157)
Y = np.linspace(8146475.91, 8160423.65, 159)
Z = np.linspace(-1082.86, -596.02, 159)
c = np.linspace(1082.86, 596.02, 159)
# X, Y = np.meshgrid(x, y)

outF = open("myOutFile_" + "gslib2sblwiz_reservoir" + ".txt", "w")
outF.write("XA" + " " + "YA" + " " + "ZA" + " " + "GA")
outF.write("\n")
for k in range(101):
    for j in range(159):
        for i in range(157):
            outF.write(str(598982.08 + (i * (11753.15 / 157))) + " " + str(8146475.91 + (j * (13947.74 / 159))) + " "
                       + str(-1082.86 + (k * (486.84 / 101))) + " " + str(converted_model[k, j, i]))
            outF.write("\n")
outF.close()

outF1 = open("myOutFile_" + "gslib2sblwiz_overburden" + ".txt", "w")
outF1.write("XA" + " " + "YA" + " " + "ZA" + " " + "GA")
outF1.write("\n")
for k in range(101):
    for j in range(159):
        for i in range(157):
            outF1.write(str(598982.08 + (i * (11753.15 / 157))) + " " + str(8146475.91 + (j * (13947.74 / 159))) + " "
                        + str(-500 + (k * (500 / 101))) + " " + str(random.randint(0, 20)))
            outF1.write("\n")

outF1.close()

outF2 = open("myOutFile_" + "gslib2sblwiz_halfspace" + ".txt", "w")
outF2.write("XA" + " " + "YA" + " " + "ZA" + " " + "GA")
outF2.write("\n")
for k in range(50):
    for j in range(159):
        for i in range(157):
            outF2.write(str(598982.08 + (i * (11753.15 / 157))) + " " + str(8146475.91 + (j * (13947.74 / 159))) + " "
                        + str(-1082.08 + (k * (300 / 50))) + " " + str(10))
            outF2.write("\n")
outF2.close()

outF3 = open("myOutFile_" + "gslib2sblwiz_complete_model" + ".txt", "w")
for k in range(200):
    for j in range(159):
        for i in range(157):
            outF3.write(str(i + 1) + " " + str(j + 1) + " " + str(k + 1) + " " + (str(x_coordinates[k, j, i]) +
                                                                                  " " + str(
                        y_coordinates[k, j, i]) + " " + str(
                        z_coordinates[k, j, i]) + " " + str(plotted_model[k, j, i])))
            outF3.write("\n")
outF.close()

LOGMIN = 0.1


########################################################################################################################
# Plot the model in 3D
# create meshgrid
# xi = np.arange(X.min(), X.max(), (X.max() - X.min()) / 100)
# yi = np.arange(Y.min(), Y.max(), (Y.max() - Y.min()) / 100)
# xi, yi = np.meshgrid(xi, yi)
# # interpolate
# zi = griddata((X, Y), Z, (xi, yi), method='linear')
# # Plot data
# fig = go.Figure(data=[
#     go.Surface(z=-ziw, x=xiw, y=yiw, opacity=1, colorscale='Greys'),
#     go.Surface(z=-zi, x=xi, y=yi, opacity=0.75, colorscale='Greens'),
#     go.Surface(z=-zi2, x=xi2, y=yi2, opacity=0.75, showscale=False, colorscale='Greens', name='Hugin Top',
#                hovertemplate="Depth: %{z:.0f}m<extra>Hugin<br>Form.</extra>"),
#     go.Surface(z=-zi3, x=xi3, y=yi3, opacity=0.75, colorscale='Greens'),
#     go.Surface(z=-zi4, x=xi4, y=yi4, opacity=0.75, colorscale='Greens')
# ])
#
# fig.update_traces(contours_z=dict(show=True, usecolormap=True,
#                                   highlightcolor="limegreen", project_z=True))
#
# zoom = 1.25
# fig.update_layout(hovermode='closest', xaxis=dict(showspikes=False),
#                   scene_camera=dict(eye=dict(x=0.5 * zoom, y=-1 * zoom, z=0.5 * zoom)))
# fig.show()


def grids_maker(filepath):
    # Get the data
    df = pd.read_csv(filepath, sep=' ')

    # Make things more legible
    xy = df[['XA', 'YA']]
    x = xy.XA
    y = xy.YA
    z = df.ZA
    g = df.GA
    reso_x = reso_y = 100
    interp = 'cubic'  # or 'nearest' or 'linear'

    # Convert the 4d-space's dimensions into grids
    grid_x, grid_y = np.mgrid[
                     x.min():x.max():1j * reso_x,
                     y.min():y.max():1j * reso_y
                     ]

    grid_z = si.griddata(
        xy, z.values,
        (grid_x, grid_y),
        method=interp
    )

    grid_g = si.griddata(
        xy, g.values,
        (grid_x, grid_y),
        method=interp
    )

    return {
        'x': grid_x,
        'y': grid_y,
        'z': grid_z,
        'g': grid_g,
    }


# Let's retrieve all files' contents
fgrids_r = dict.fromkeys([
    'myOutFile_gslib2sblwiz_reservoir.txt'
])
g_r_mins = []
g_r_maxs = []

fgrids_o = dict.fromkeys([
    'myOutFile_gslib2sblwiz_overburden.txt'
])
g_o_mins = []
g_o_maxs = []

fgrids_h = dict.fromkeys([
    'myOutFile_gslib2sblwiz_halfspace.txt'
])
g_h_mins = []
g_h_maxs = []

for fpath in fgrids_r.keys():
    fgrids_r[fpath] = grids_r = grids_maker('myOutFile_gslib2sblwiz_reservoir.txt')
    g_r_mins.append(grids_r['g'].min())
    g_r_maxs.append(grids_r['g'].max())

result = np.zeros((10, 100, 100))
for i in range(10):
    for fpath in fgrids_o.keys():
        fgrids_o[fpath] = grids_o = grids_maker('myOutFile_gslib2sblwiz_overburden.txt')
        g_o_mins.append(grids_o['g'].min())
        g_o_maxs.append(grids_o['g'].max())
        result[i] = grids_o['z']

for fpath in fgrids_h.keys():
    fgrids_h[fpath] = grids_h = grids_maker('myOutFile_gslib2sblwiz_halfspace.txt')
    g_h_mins.append(grids_h['g'].min())
    g_h_maxs.append(grids_h['g'].max())

# Create the 4th color-rendered dimension
scam_r = plt.cm.ScalarMappable(
    norm=cm.colors.Normalize(min(g_r_mins), max(g_r_maxs)),
    cmap='jet'  # see https://matplotlib.org/examples/color/colormaps_reference.html
)

# Create the 4th color-rendered dimension
scam_o = plt.cm.ScalarMappable(
    norm=cm.colors.Normalize(min(g_o_mins), max(g_o_maxs)),
    cmap='jet'  # see https://matplotlib.org/examples/color/colormaps_reference.html
)

# Create the 4th color-rendered dimension
scam_h = plt.cm.ScalarMappable(
    norm=cm.colors.Normalize(min(g_h_mins), max(g_h_maxs)),
    cmap='jet'  # see https://matplotlib.org/examples/color/colormaps_reference.html
)

# Make the plot
fig = plt.figure()
ax = fig.gca(projection='3d')

# # Reservoir
# for grids_r in fgrids_r.values():
#     scam_r.set_array([])
#     surf_r = ax.plot_surface(
#         grids_r['x'], grids_r['y'], grids_r['z'],
#         facecolors=scam_r.to_rgba(grids_r['g']),
#         antialiased=True,
#         rstride=1, cstride=1, alpha=None

# )
# Overburden
for i in range(10):
    for grids_o in fgrids_o.values():
        grids_o['z'] = result[i]
        scam_o.set_array([])
        surf_o = ax.plot_surface(
            grids_o['x'], grids_o['y'], grids_o['z'],
            facecolors=scam_o.to_rgba(grids_o['g']),
            antialiased=True,
            rstride=1, cstride=1, alpha=None

        )
        #     # Halfspace
        # for grids_h in fgrids_h.values():
        #     scam_h.set_array([])
        #     surf_h = ax.plot_surface(
        #         grids_h['x'], grids_h['y'], grids_h['z'],
        #         facecolors=scam_h.to_rgba(grids_h['g']),
        #         antialiased=True,
        #         rstride=1, cstride=1, alpha=None

        # )
        # Customize the z axis.
        ax.set_zlim(-10000, 0)
        fig.colorbar(surf_o, shrink=0.5, aspect=5)
    plt.show()

# converted_model = converted_model.flatten()
#
# # triangles = mtri.Triangulation(X, Y).triangles
#
#
# # choice_calcuation_colors = 1
# # if choice_calcuation_colors == 1:  # Mean of the "c" values of the 3 pt of the triangle
# #     colors = np.mean([c[triangles[:, 0]], c[triangles[:, 1]], c[triangles[:, 2]]], axis=0)
# # elif choice_calcuation_colors == 2:  # Mediane of the "c" values of the 3 pt of the triangle
# #     colors = np.median([c[triangles[:, 0]], c[triangles[:, 1]], c[triangles[:, 2]]], axis=0)
# # elif choice_calcuation_colors == 3:  # Max of the "c" values of the 3 pt of the triangle
# #     colors = np.max([c[triangles[:, 0]], c[triangles[:, 1]], c[triangles[:, 2]]], axis=0)
#
# # Displays the 4D graphic.
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# img = ax.scatter(X, Y, Z, c=converted_model, cmap=plt.hot())
# fig.colorbar(img)
# plt.show()
#
# # Add a color bar with a title to explain which variable is represented by the color.
# # cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
# # cbar.ax.get_yaxis().labelpad = 15
# # cbar.ax.set_ylabel(list_name_variables[index_c], rotation=270)
# #
# # # Add titles to the axes and a title in the figure.
# # ax.set_xlabel(list_name_variables[index_x])
# # ax.set_ylabel(list_name_variables[index_y])
# # ax.set_zlabel(list_name_variables[index_Z])
# # plt.title('%s in function of %s, %s and %s' % (
# #     list_name_variables[index_c], list_name_variables[index_x], list_name_variables[index_y],
# #     list_name_variables[index_z]))

# i = 0
# while i < 3:
M = 0
N = 0
for i in range(200):
    M += z_coordinates[i, :, :]
    N += plotted_model[i, :, :]
M = M / 200
N = N / 200
fig = go.Figure(data=[
    go.Surface(z=z_coordinates[1, :, :], surfacecolor=plotted_model[1, :, :], cmin=0, cmax=35, colorscale='Jet'),
    go.Surface(z=z_coordinates[2, :, :], surfacecolor=plotted_model[2, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[3, :, :], surfacecolor=plotted_model[3, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[4, :, :], surfacecolor=plotted_model[4, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[5, :, :], surfacecolor=plotted_model[5, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[6, :, :], surfacecolor=plotted_model[6, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[7, :, :], surfacecolor=plotted_model[7, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[8, :, :], surfacecolor=plotted_model[8, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[9, :, :], surfacecolor=plotted_model[9, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[10, :, :], surfacecolor=plotted_model[10, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[11:, :], surfacecolor=plotted_model[11:, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[12:, :], surfacecolor=plotted_model[12:, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[13, :, :], surfacecolor=plotted_model[13, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[14, :, :], surfacecolor=plotted_model[14, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[15, :, :], surfacecolor=plotted_model[15, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[16, :, :], surfacecolor=plotted_model[16, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[17, :, :], surfacecolor=plotted_model[17, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[18, :, :], surfacecolor=plotted_model[18, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[19, :, :], surfacecolor=plotted_model[19, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[20, :, :], surfacecolor=plotted_model[20, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[21, :, :], surfacecolor=plotted_model[21, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[22, :, :], surfacecolor=plotted_model[22, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[23, :, :], surfacecolor=plotted_model[23, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[24, :, :], surfacecolor=plotted_model[24, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[25, :, :], surfacecolor=plotted_model[25, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[26, :, :], surfacecolor=plotted_model[26, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[27, :, :], surfacecolor=plotted_model[27, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[28, :, :], surfacecolor=plotted_model[28, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[29, :, :], surfacecolor=plotted_model[29, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[30, :, :], surfacecolor=plotted_model[30, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[31, :, :], surfacecolor=plotted_model[31, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[32, :, :], surfacecolor=plotted_model[32, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[33, :, :], surfacecolor=plotted_model[33, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[34, :, :], surfacecolor=plotted_model[34, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[35, :, :], surfacecolor=plotted_model[35, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[36, :, :], surfacecolor=plotted_model[36, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[37, :, :], surfacecolor=plotted_model[37, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[38, :, :], surfacecolor=plotted_model[38, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[39, :, :], surfacecolor=plotted_model[39, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[40, :, :], surfacecolor=plotted_model[40, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[41, :, :], surfacecolor=plotted_model[41, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[42, :, :], surfacecolor=plotted_model[42, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[43, :, :], surfacecolor=plotted_model[43, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[44, :, :], surfacecolor=plotted_model[44, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[45, :, :], surfacecolor=plotted_model[45, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[46, :, :], surfacecolor=plotted_model[46, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[47, :, :], surfacecolor=plotted_model[47, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[48, :, :], surfacecolor=plotted_model[48, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[49, :, :], surfacecolor=plotted_model[49, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[50, :, :], surfacecolor=plotted_model[50, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[51, :, :], surfacecolor=plotted_model[51, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[52, :, :], surfacecolor=plotted_model[52, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[53, :, :], surfacecolor=plotted_model[53, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[54, :, :], surfacecolor=plotted_model[54, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[55, :, :], surfacecolor=plotted_model[55, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[56, :, :], surfacecolor=plotted_model[56, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[57, :, :], surfacecolor=plotted_model[57, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[58, :, :], surfacecolor=plotted_model[58, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[59, :, :], surfacecolor=plotted_model[59, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[60, :, :], surfacecolor=plotted_model[60, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[61, :, :], surfacecolor=plotted_model[61, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[62, :, :], surfacecolor=plotted_model[62, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[63, :, :], surfacecolor=plotted_model[63, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[64, :, :], surfacecolor=plotted_model[64, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[65, :, :], surfacecolor=plotted_model[65, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[66, :, :], surfacecolor=plotted_model[66, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[67, :, :], surfacecolor=plotted_model[67, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[68, :, :], surfacecolor=plotted_model[68, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[69, :, :], surfacecolor=plotted_model[69, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[70, :, :], surfacecolor=plotted_model[70, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[71, :, :], surfacecolor=plotted_model[71, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[72, :, :], surfacecolor=plotted_model[72, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[73, :, :], surfacecolor=plotted_model[73, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[74, :, :], surfacecolor=plotted_model[74, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[75, :, :], surfacecolor=plotted_model[75, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[76, :, :], surfacecolor=plotted_model[76, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[77, :, :], surfacecolor=plotted_model[77, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[78, :, :], surfacecolor=plotted_model[78, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[79, :, :], surfacecolor=plotted_model[79, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[80, :, :], surfacecolor=plotted_model[80, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[81, :, :], surfacecolor=plotted_model[81, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[82, :, :], surfacecolor=plotted_model[82, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[83, :, :], surfacecolor=plotted_model[83, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[84, :, :], surfacecolor=plotted_model[84, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[85, :, :], surfacecolor=plotted_model[85, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[86, :, :], surfacecolor=plotted_model[86, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[87, :, :], surfacecolor=plotted_model[87, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[88, :, :], surfacecolor=plotted_model[88, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[89, :, :], surfacecolor=plotted_model[89, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[90, :, :], surfacecolor=plotted_model[90, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[91, :, :], surfacecolor=plotted_model[91, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[92, :, :], surfacecolor=plotted_model[92, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[93, :, :], surfacecolor=plotted_model[93, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[94, :, :], surfacecolor=plotted_model[94, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[95, :, :], surfacecolor=plotted_model[95, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[96, :, :], surfacecolor=plotted_model[96, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[97, :, :], surfacecolor=plotted_model[97, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[98, :, :], surfacecolor=plotted_model[98, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[99, :, :], surfacecolor=plotted_model[99, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[100, :, :], surfacecolor=plotted_model[100, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[101, :, :], surfacecolor=plotted_model[101, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[102, :, :], surfacecolor=plotted_model[102, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[103, :, :], surfacecolor=plotted_model[103, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[104, :, :], surfacecolor=plotted_model[104, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[105, :, :], surfacecolor=plotted_model[105, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[106:, :], surfacecolor=plotted_model[106, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[107:, :], surfacecolor=plotted_model[107, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[108:, :], surfacecolor=plotted_model[108, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[109, :, :], surfacecolor=plotted_model[109, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[110, :, :], surfacecolor=plotted_model[110, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[111, :, :], surfacecolor=plotted_model[111, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[112, :, :], surfacecolor=plotted_model[112, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[113, :, :], surfacecolor=plotted_model[113, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[114, :, :], surfacecolor=plotted_model[114, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[115, :, :], surfacecolor=plotted_model[115, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[116, :, :], surfacecolor=plotted_model[116, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[117, :, :], surfacecolor=plotted_model[117, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[118, :, :], surfacecolor=plotted_model[118, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[119, :, :], surfacecolor=plotted_model[119, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[120, :, :], surfacecolor=plotted_model[120, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[121, :, :], surfacecolor=plotted_model[121, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[122, :, :], surfacecolor=plotted_model[122, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[123, :, :], surfacecolor=plotted_model[123, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[124, :, :], surfacecolor=plotted_model[124, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[125, :, :], surfacecolor=plotted_model[125, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[126, :, :], surfacecolor=plotted_model[126, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[127, :, :], surfacecolor=plotted_model[127, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[128, :, :], surfacecolor=plotted_model[128, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[129, :, :], surfacecolor=plotted_model[129, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[130, :, :], surfacecolor=plotted_model[130, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[131, :, :], surfacecolor=plotted_model[131, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[132, :, :], surfacecolor=plotted_model[132, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[133, :, :], surfacecolor=plotted_model[133, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[134, :, :], surfacecolor=plotted_model[134, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[135, :, :], surfacecolor=plotted_model[135, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[136, :, :], surfacecolor=plotted_model[136, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[137, :, :], surfacecolor=plotted_model[137, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[138, :, :], surfacecolor=plotted_model[138, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[139, :, :], surfacecolor=plotted_model[139, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[140, :, :], surfacecolor=plotted_model[140, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[141, :, :], surfacecolor=plotted_model[141, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[142, :, :], surfacecolor=plotted_model[142, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[143, :, :], surfacecolor=plotted_model[143, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[144, :, :], surfacecolor=plotted_model[144, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[145, :, :], surfacecolor=plotted_model[145, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[146, :, :], surfacecolor=plotted_model[146, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[147, :, :], surfacecolor=plotted_model[147, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[148, :, :], surfacecolor=plotted_model[148, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[149, :, :], surfacecolor=plotted_model[149, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[150, :, :], surfacecolor=plotted_model[150, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[151, :, :], surfacecolor=plotted_model[151, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[152, :, :], surfacecolor=plotted_model[152, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[153, :, :], surfacecolor=plotted_model[153, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[154, :, :], surfacecolor=plotted_model[154, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[155, :, :], surfacecolor=plotted_model[155, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[156, :, :], surfacecolor=plotted_model[156, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[157, :, :], surfacecolor=plotted_model[157, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[158, :, :], surfacecolor=plotted_model[158, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[159, :, :], surfacecolor=plotted_model[159, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[160, :, :], surfacecolor=plotted_model[160, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[161, :, :], surfacecolor=plotted_model[161, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[162, :, :], surfacecolor=plotted_model[162, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[163, :, :], surfacecolor=plotted_model[163, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[164, :, :], surfacecolor=plotted_model[164, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[165, :, :], surfacecolor=plotted_model[165, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[166, :, :], surfacecolor=plotted_model[166, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[167, :, :], surfacecolor=plotted_model[167, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[168, :, :], surfacecolor=plotted_model[168, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[169, :, :], surfacecolor=plotted_model[169, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[170, :, :], surfacecolor=plotted_model[170, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[171, :, :], surfacecolor=plotted_model[171, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[172, :, :], surfacecolor=plotted_model[172, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[173, :, :], surfacecolor=plotted_model[173, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[174, :, :], surfacecolor=plotted_model[174, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[175, :, :], surfacecolor=plotted_model[175, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[176, :, :], surfacecolor=plotted_model[176, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[177, :, :], surfacecolor=plotted_model[177, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[178, :, :], surfacecolor=plotted_model[178, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[179, :, :], surfacecolor=plotted_model[179, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[180, :, :], surfacecolor=plotted_model[180, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[181, :, :], surfacecolor=plotted_model[181, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[182, :, :], surfacecolor=plotted_model[182, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[183, :, :], surfacecolor=plotted_model[183, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[184, :, :], surfacecolor=plotted_model[184, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[185, :, :], surfacecolor=plotted_model[185, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[186, :, :], surfacecolor=plotted_model[186, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[187, :, :], surfacecolor=plotted_model[187, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[188, :, :], surfacecolor=plotted_model[188, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[189, :, :], surfacecolor=plotted_model[189, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[190, :, :], surfacecolor=plotted_model[190, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[191, :, :], surfacecolor=plotted_model[191, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[192, :, :], surfacecolor=plotted_model[192, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[193, :, :], surfacecolor=plotted_model[193, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[194, :, :], surfacecolor=plotted_model[194, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[195, :, :], surfacecolor=plotted_model[195, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[196, :, :], surfacecolor=plotted_model[196, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[197, :, :], surfacecolor=plotted_model[197, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[198, :, :], surfacecolor=plotted_model[198, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[199, :, :], surfacecolor=plotted_model[199, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[200, :, :], surfacecolor=plotted_model[200, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[201, :, :], surfacecolor=plotted_model[201, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[202, :, :], surfacecolor=plotted_model[202, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[203, :, :], surfacecolor=plotted_model[203, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[204, :, :], surfacecolor=plotted_model[204, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[205, :, :], surfacecolor=plotted_model[205, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[206, :, :], surfacecolor=plotted_model[206, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[207, :, :], surfacecolor=plotted_model[207, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[208, :, :], surfacecolor=plotted_model[208, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[209, :, :], surfacecolor=plotted_model[209, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[210, :, :], surfacecolor=plotted_model[210, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[211, :, :], surfacecolor=plotted_model[211, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[212, :, :], surfacecolor=plotted_model[212, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[213, :, :], surfacecolor=plotted_model[213, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[214, :, :], surfacecolor=plotted_model[214, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[215, :, :], surfacecolor=plotted_model[215, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[216, :, :], surfacecolor=plotted_model[216, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[217, :, :], surfacecolor=plotted_model[217, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[218, :, :], surfacecolor=plotted_model[218, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[219, :, :], surfacecolor=plotted_model[219, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[220, :, :], surfacecolor=plotted_model[220, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[221, :, :], surfacecolor=plotted_model[221, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[222, :, :], surfacecolor=plotted_model[222, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[223, :, :], surfacecolor=plotted_model[223, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[224, :, :], surfacecolor=plotted_model[224, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[225, :, :], surfacecolor=plotted_model[225, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[226, :, :], surfacecolor=plotted_model[226, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[227, :, :], surfacecolor=plotted_model[227, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[228, :, :], surfacecolor=plotted_model[228, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[229, :, :], surfacecolor=plotted_model[229, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[230, :, :], surfacecolor=plotted_model[230, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[231, :, :], surfacecolor=plotted_model[231, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[232, :, :], surfacecolor=plotted_model[232, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[233, :, :], surfacecolor=plotted_model[233, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[234, :, :], surfacecolor=plotted_model[234, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[235, :, :], surfacecolor=plotted_model[235, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[236, :, :], surfacecolor=plotted_model[236, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[237, :, :], surfacecolor=plotted_model[237, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[238, :, :], surfacecolor=plotted_model[238, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[239, :, :], surfacecolor=plotted_model[239, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[240, :, :], surfacecolor=plotted_model[240, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[241, :, :], surfacecolor=plotted_model[241, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[242, :, :], surfacecolor=plotted_model[242, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[243, :, :], surfacecolor=plotted_model[243, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[244, :, :], surfacecolor=plotted_model[244, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[245, :, :], surfacecolor=plotted_model[245, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[246, :, :], surfacecolor=plotted_model[246, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[247, :, :], surfacecolor=plotted_model[247, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[248, :, :], surfacecolor=plotted_model[248, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[249, :, :], surfacecolor=plotted_model[249, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[250, :, :], surfacecolor=plotted_model[250, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[251, :, :], surfacecolor=plotted_model[251, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[252, :, :], surfacecolor=plotted_model[252, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[253, :, :], surfacecolor=plotted_model[253, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[254, :, :], surfacecolor=plotted_model[254, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[255, :, :], surfacecolor=plotted_model[255, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[256, :, :], surfacecolor=plotted_model[256, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[257, :, :], surfacecolor=plotted_model[257, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[258, :, :], surfacecolor=plotted_model[258, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[259, :, :], surfacecolor=plotted_model[259, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[260, :, :], surfacecolor=plotted_model[260, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[261, :, :], surfacecolor=plotted_model[261, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[262, :, :], surfacecolor=plotted_model[262, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[263, :, :], surfacecolor=plotted_model[263, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[264, :, :], surfacecolor=plotted_model[264, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[265, :, :], surfacecolor=plotted_model[265, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[266, :, :], surfacecolor=plotted_model[266, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[267, :, :], surfacecolor=plotted_model[267, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[268, :, :], surfacecolor=plotted_model[268, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[269, :, :], surfacecolor=plotted_model[269, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[270, :, :], surfacecolor=plotted_model[270, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[271, :, :], surfacecolor=plotted_model[271, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[272, :, :], surfacecolor=plotted_model[272, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[273, :, :], surfacecolor=plotted_model[273, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[274, :, :], surfacecolor=plotted_model[274, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[275, :, :], surfacecolor=plotted_model[275, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[276, :, :], surfacecolor=plotted_model[276, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[277, :, :], surfacecolor=plotted_model[277, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[278, :, :], surfacecolor=plotted_model[278, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[279, :, :], surfacecolor=plotted_model[279, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[280, :, :], surfacecolor=plotted_model[280, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[281, :, :], surfacecolor=plotted_model[281, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[282, :, :], surfacecolor=plotted_model[282, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[283, :, :], surfacecolor=plotted_model[283, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[284, :, :], surfacecolor=plotted_model[284, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[285, :, :], surfacecolor=plotted_model[285, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[286, :, :], surfacecolor=plotted_model[286, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[287, :, :], surfacecolor=plotted_model[287, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[288, :, :], surfacecolor=plotted_model[288, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[289, :, :], surfacecolor=plotted_model[289, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[290, :, :], surfacecolor=plotted_model[290, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[291, :, :], surfacecolor=plotted_model[291, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[292, :, :], surfacecolor=plotted_model[292, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[293, :, :], surfacecolor=plotted_model[293, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[294, :, :], surfacecolor=plotted_model[294, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[295, :, :], surfacecolor=plotted_model[295, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[296, :, :], surfacecolor=plotted_model[296, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[297, :, :], surfacecolor=plotted_model[297, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[298, :, :], surfacecolor=plotted_model[298, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[299, :, :], surfacecolor=plotted_model[299, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[300, :, :], surfacecolor=plotted_model[300, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[301, :, :], surfacecolor=plotted_model[301, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[302, :, :], surfacecolor=plotted_model[302, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[303, :, :], surfacecolor=plotted_model[303, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[304, :, :], surfacecolor=plotted_model[304, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[305, :, :], surfacecolor=plotted_model[305, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[306, :, :], surfacecolor=plotted_model[306, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[307, :, :], surfacecolor=plotted_model[307, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[308, :, :], surfacecolor=plotted_model[308, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[309, :, :], surfacecolor=plotted_model[309, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[311, :, :], surfacecolor=plotted_model[311, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[312, :, :], surfacecolor=plotted_model[312, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[312, :, :], surfacecolor=plotted_model[312, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[314, :, :], surfacecolor=plotted_model[314, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[315, :, :], surfacecolor=plotted_model[315, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[316, :, :], surfacecolor=plotted_model[316, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[317, :, :], surfacecolor=plotted_model[317, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[318, :, :], surfacecolor=plotted_model[318, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[319, :, :], surfacecolor=plotted_model[319, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[320, :, :], surfacecolor=plotted_model[320, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[321, :, :], surfacecolor=plotted_model[321, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[322, :, :], surfacecolor=plotted_model[322, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[323, :, :], surfacecolor=plotted_model[323, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[324, :, :], surfacecolor=plotted_model[324, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[325, :, :], surfacecolor=plotted_model[325, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[326, :, :], surfacecolor=plotted_model[326, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[327, :, :], surfacecolor=plotted_model[327, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[328, :, :], surfacecolor=plotted_model[328, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[329, :, :], surfacecolor=plotted_model[329, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[330, :, :], surfacecolor=plotted_model[330, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[331, :, :], surfacecolor=plotted_model[331, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[332, :, :], surfacecolor=plotted_model[332, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[333, :, :], surfacecolor=plotted_model[333, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[334, :, :], surfacecolor=plotted_model[334, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[335, :, :], surfacecolor=plotted_model[335, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[336, :, :], surfacecolor=plotted_model[336, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[337, :, :], surfacecolor=plotted_model[337, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[338, :, :], surfacecolor=plotted_model[338, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[339, :, :], surfacecolor=plotted_model[339, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[340, :, :], surfacecolor=plotted_model[340, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[341, :, :], surfacecolor=plotted_model[341, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[342, :, :], surfacecolor=plotted_model[342, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[343, :, :], surfacecolor=plotted_model[343, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[344, :, :], surfacecolor=plotted_model[344, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[345, :, :], surfacecolor=plotted_model[345, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[346, :, :], surfacecolor=plotted_model[346, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[347, :, :], surfacecolor=plotted_model[347, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[348, :, :], surfacecolor=plotted_model[348, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[349, :, :], surfacecolor=plotted_model[349, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[350, :, :], surfacecolor=plotted_model[350, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[351, :, :], surfacecolor=plotted_model[351, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[352, :, :], surfacecolor=plotted_model[352, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[353, :, :], surfacecolor=plotted_model[353, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[354, :, :], surfacecolor=plotted_model[354, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[355, :, :], surfacecolor=plotted_model[355, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[356, :, :], surfacecolor=plotted_model[356, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[357, :, :], surfacecolor=plotted_model[357, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[358, :, :], surfacecolor=plotted_model[358, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[359, :, :], surfacecolor=plotted_model[359, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[360, :, :], surfacecolor=plotted_model[360, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[361, :, :], surfacecolor=plotted_model[361, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[362, :, :], surfacecolor=plotted_model[362, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[363, :, :], surfacecolor=plotted_model[363, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[364, :, :], surfacecolor=plotted_model[364, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[365, :, :], surfacecolor=plotted_model[365, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[366, :, :], surfacecolor=plotted_model[366, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[367, :, :], surfacecolor=plotted_model[367, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[368, :, :], surfacecolor=plotted_model[368, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[369, :, :], surfacecolor=plotted_model[369, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[370, :, :], surfacecolor=plotted_model[370, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[371, :, :], surfacecolor=plotted_model[371, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[372, :, :], surfacecolor=plotted_model[372, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[373, :, :], surfacecolor=plotted_model[373, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[374, :, :], surfacecolor=plotted_model[374, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[375, :, :], surfacecolor=plotted_model[375, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[376, :, :], surfacecolor=plotted_model[376, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[377, :, :], surfacecolor=plotted_model[377, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[378, :, :], surfacecolor=plotted_model[378, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[379, :, :], surfacecolor=plotted_model[379, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[380, :, :], surfacecolor=plotted_model[380, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[381, :, :], surfacecolor=plotted_model[381, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[382, :, :], surfacecolor=plotted_model[382, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[383, :, :], surfacecolor=plotted_model[383, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[384, :, :], surfacecolor=plotted_model[384, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[385, :, :], surfacecolor=plotted_model[385, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[386, :, :], surfacecolor=plotted_model[386, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[387, :, :], surfacecolor=plotted_model[387, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[388, :, :], surfacecolor=plotted_model[388, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[389, :, :], surfacecolor=plotted_model[389, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[390, :, :], surfacecolor=plotted_model[390, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[391, :, :], surfacecolor=plotted_model[391, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[392, :, :], surfacecolor=plotted_model[392, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[393, :, :], surfacecolor=plotted_model[393, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[394, :, :], surfacecolor=plotted_model[394, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[395, :, :], surfacecolor=plotted_model[395, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[396, :, :], surfacecolor=plotted_model[396, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[397, :, :], surfacecolor=plotted_model[397, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[398, :, :], surfacecolor=plotted_model[398, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[399, :, :], surfacecolor=plotted_model[399, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[400, :, :], surfacecolor=plotted_model[400, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[401, :, :], surfacecolor=plotted_model[401, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[402, :, :], surfacecolor=plotted_model[402, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[403, :, :], surfacecolor=plotted_model[403, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[404, :, :], surfacecolor=plotted_model[404, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[405, :, :], surfacecolor=plotted_model[405, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[406, :, :], surfacecolor=plotted_model[406, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[407, :, :], surfacecolor=plotted_model[407, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[408, :, :], surfacecolor=plotted_model[408, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[409, :, :], surfacecolor=plotted_model[409, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[410, :, :], surfacecolor=plotted_model[410, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[411, :, :], surfacecolor=plotted_model[411, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[412, :, :], surfacecolor=plotted_model[412, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[413, :, :], surfacecolor=plotted_model[413, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[414, :, :], surfacecolor=plotted_model[414, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[415, :, :], surfacecolor=plotted_model[415, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[416, :, :], surfacecolor=plotted_model[416, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[417, :, :], surfacecolor=plotted_model[417, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[418, :, :], surfacecolor=plotted_model[418, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[419, :, :], surfacecolor=plotted_model[419, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[420, :, :], surfacecolor=plotted_model[420, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[421, :, :], surfacecolor=plotted_model[421, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[422, :, :], surfacecolor=plotted_model[422, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[423, :, :], surfacecolor=plotted_model[423, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[424, :, :], surfacecolor=plotted_model[424, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[425, :, :], surfacecolor=plotted_model[425, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[426, :, :], surfacecolor=plotted_model[426, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[427, :, :], surfacecolor=plotted_model[427, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[428, :, :], surfacecolor=plotted_model[428, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[429, :, :], surfacecolor=plotted_model[429, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[430, :, :], surfacecolor=plotted_model[430, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[431, :, :], surfacecolor=plotted_model[431, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[432, :, :], surfacecolor=plotted_model[432, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[433, :, :], surfacecolor=plotted_model[433, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[434, :, :], surfacecolor=plotted_model[434, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[435, :, :], surfacecolor=plotted_model[435, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[436, :, :], surfacecolor=plotted_model[436, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[437, :, :], surfacecolor=plotted_model[437, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[438, :, :], surfacecolor=plotted_model[438, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[439, :, :], surfacecolor=plotted_model[439, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[440, :, :], surfacecolor=plotted_model[440, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[441, :, :], surfacecolor=plotted_model[441, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[442, :, :], surfacecolor=plotted_model[442, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[443, :, :], surfacecolor=plotted_model[443, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[444, :, :], surfacecolor=plotted_model[444, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[445, :, :], surfacecolor=plotted_model[445, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[446, :, :], surfacecolor=plotted_model[446, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[447, :, :], surfacecolor=plotted_model[447, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[448, :, :], surfacecolor=plotted_model[448, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[449, :, :], surfacecolor=plotted_model[449, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[450, :, :], surfacecolor=plotted_model[450, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[451, :, :], surfacecolor=plotted_model[451, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[452, :, :], surfacecolor=plotted_model[452, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[453, :, :], surfacecolor=plotted_model[453, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[454, :, :], surfacecolor=plotted_model[454, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[455, :, :], surfacecolor=plotted_model[455, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[456, :, :], surfacecolor=plotted_model[456, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[457, :, :], surfacecolor=plotted_model[457, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[458, :, :], surfacecolor=plotted_model[458, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[459, :, :], surfacecolor=plotted_model[459, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[460, :, :], surfacecolor=plotted_model[460, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[461, :, :], surfacecolor=plotted_model[461, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[462, :, :], surfacecolor=plotted_model[462, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[463, :, :], surfacecolor=plotted_model[463, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[464, :, :], surfacecolor=plotted_model[464, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[465, :, :], surfacecolor=plotted_model[465, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[466, :, :], surfacecolor=plotted_model[466, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[467, :, :], surfacecolor=plotted_model[467, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[468, :, :], surfacecolor=plotted_model[468, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[469, :, :], surfacecolor=plotted_model[469, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[470, :, :], surfacecolor=plotted_model[470, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[471, :, :], surfacecolor=plotted_model[471, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[472, :, :], surfacecolor=plotted_model[472, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[473, :, :], surfacecolor=plotted_model[473, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[474, :, :], surfacecolor=plotted_model[474, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[475, :, :], surfacecolor=plotted_model[475, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[476, :, :], surfacecolor=plotted_model[476, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[477, :, :], surfacecolor=plotted_model[477, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[478, :, :], surfacecolor=plotted_model[478, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[479, :, :], surfacecolor=plotted_model[479, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[480, :, :], surfacecolor=plotted_model[480, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[481, :, :], surfacecolor=plotted_model[481, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[482, :, :], surfacecolor=plotted_model[482, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[483, :, :], surfacecolor=plotted_model[483, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[484, :, :], surfacecolor=plotted_model[484, :, :], cmin=0, cmax=35, colorscale='jet'),
    go.Surface(z=z_coordinates[485, :, :], surfacecolor=plotted_model[485, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
    # go.Surface(z=z_coordinates[310, :, :], surfacecolor=plotted_model[310, :, :], cmin=0, cmax=35, colorscale='jet'),
])
fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))
# i += 1
zoom = 1.25
fig.update_layout(hovermode='closest', xaxis=dict(showspikes=False),
                  scene_camera=dict(eye=dict(x=0.5 * zoom, y=-1 * zoom, z=0.5 * zoom)))
fig.show()

# if __name__ == '__main__':
#     app.run_server(debug=True, use_reloader=False, port=8051)
