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
            k = int(arr[2]) - 1
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
    model = np.zeros((101, 159, 157))
    x = np.zeros((101, 159, 157))
    y = np.zeros((101, 159, 157))
    z = np.zeros((101, 159, 157))

    with open(file, "r") as file:
        lines = file.readlines()
        file.close()

    # del lines[0]

    for line in lines:
        arr = line.split()
        if line != "\n":
            line_count += 1
        index = 0

        for s in range(len(arr)):

            i = int(arr[0]) - 1
            index += 1
            j = int(arr[1]) - 1
            index += 1
            k = int(arr[2]) - 1
            index += 1
            model[k, j, i] = float(arr[6])
            x[k, j, i] = float(arr[3])
            y[k, j, i] = float(arr[4])
            z[k, j, i] = float(arr[5])
            if index == 3:
                index = 0
                continue

    file.close()

    return model, x, y, z


########################################################################################################################


########################################################################################################################
resis_target = 'moon.txt'
resis_target_edited = 'moon_edited.txt'
converted_model = petrel2sblwiz(resis_target, resis_target_edited)
plotted_model, x_coordinates, y_coordinates, z_coordinates = multiple3dsurfaceplots(resis_target_edited)
# handles = [Rectangle((0, 0), 1, 1, color=Color, ec="k")]
simple_model, x_simple, y_simple, z_simple = simplemodelmaker()

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

LOGMIN = 0.1


########################################################################################################################


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
for i in range(101):
    M += z_coordinates[i, :, :]
    N += plotted_model[i, :, :]
M = M / 101
N = N / 101
fig = go.Figure(data=[
    go.Surface(z=M, surfacecolor=N, cmin=0, cmax=35, colorscale='Jet'),
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
