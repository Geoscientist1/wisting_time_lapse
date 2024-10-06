import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import pandas as pd

# sat = np.loadtxt("sat57_coarse.txt")
# points = np.array([[0, 0, 0], [1, 0, 0], [1, 0.5, 0], [0, 0.5, 0]])
# pdata = pyvista.PolyData(points)
# pdata['data'] = np.arange(len(points))
#
#
# sphere = pyvista.Sphere(radius=0.02, phi_resolution=10, theta_resolution=10)
# pc = pdata.glyph(scale=False, geom=sphere, orient=False)
# pc.plot(cmap='Reds')


sat27 = np.loadtxt('sat27_coarse.txt', delimiter=' ')
sat57 = np.loadtxt('sat57_coarse.txt', delimiter=' ')
ref_model = np.loadtxt('reference_model.txt', delimiter=' ')
dt = np.loadtxt('start_model_117.txt', delimiter=' ')
dswat_observed = np.subtract(sat57, sat27)
dswat_recovered = np.subtract(dt, sat27)

sat = "C:/Users/moham/OneDrive - NTNU/PhD/The wisting Project/Modeling/Gavimetri/sat27_coarse.txt"
# dt = pd.read_csv(sat, header=None, delim_whitespace=True, names=["a", "b", "c", "d"])
# dt = pd.read_csv(sat, header=None, delim_whitespace=True, names=["a"])
# vertices = dt.values[:, 0:3]
# cell_values = dt.values[:, 3]
# cell_values = dt.values[:]
# cell_values = dt
cell_values = dt
# for i in range(len(cell_values)):
#     if cell_values[i] > 0.45:
#         cell_values[i] = 0
# dt = np.reshape(dt, (5, 159, 157))
# dt0 = np.reshape(dt0, (5, 159, 157))
# cell_values = dt0[0]
xrng = np.linspace(13213.48, 0, 157)
yrng = np.linspace(0, 14173.40, 159)
zrng = np.linspace(-596.02, -1082.86, 5)
grid = pv.RectilinearGrid(xrng, yrng, zrng)
########################################################################################################################
plotter = pv.Plotter()
sargs = dict(height=0.25, vertical=True, position_x=0.80, position_y=0.05, title_font_size=18,
             title="Change in Water Saturation", label_font_size=20, shadow=True, n_labels=3, italic=False,
             fmt="%.1f", font_family="arial")
orig_map = plt.cm.get_cmap('jet')
reversed_map = orig_map.reversed()
plotter.background_color = "white"
plotter.add_axes(line_width=5)
plotter.hide_axes()
plotter.set_scale(zscale=2)
plotter.camera.view_angle = 60.0
plotter.add_bounding_box(line_width=15, color='black', outline=True)
actor = plotter.add_mesh(grid, show_edges=False, edge_opacity=0.1, edge_color="black", cmap=orig_map, clim=(0.0, 0.5),
                         scalars=cell_values, scalar_bar_args=sargs, reset_camera=True, render=True)
plotter.show()
########################################################################################################################
# cpos = plotter.camera_position
cpos = [(27873.672636091094, 27819.16571035057, 22089.181844846265),
        (6606.24, 7086.2, -1081.86),
        (-0.4359837926675337, -0.4339254547130702, 0.7884331501676527)]
print(cpos)
plotter = pv.Plotter(off_screen=True)
plotter.camera_position = cpos
sargs = dict(height=0.25, vertical=True, position_x=1.0, position_y=0.05, title_font_size=18,
             title="Change in Water Saturation", label_font_size=20, shadow=True, n_labels=3, italic=False,
             fmt="%.1f", font_family="arial")
orig_map = plt.cm.get_cmap('jet')
reversed_map = orig_map.reversed()
plotter.background_color = "white"
plotter.hide_axes()
plotter.set_scale(zscale=2)
# pl.camera_position = 'yz'
# plotter.camera.view_angle = 60.0
plotter.add_bounding_box(line_width=27, color='black')
plotter.add_mesh(grid, show_edges=False, edge_opacity=0.1, edge_color="black", cmap=orig_map,
                 scalars=cell_values, scalar_bar_args=sargs, show_scalar_bar=False, reset_camera=True)

# plotter.show_bounds(
#     grid='front',
#     location='outer',
#     ticks='both',
#     show_xlabels=True,
#     show_ylabels=True,
#     show_zlabels=True,
#     xtitle='Easting',
#     ytitle='Northing',
#     ztitle='Depth',
#     n_zlabels=2,
#     font_size=17,
# )

plotter.show(screenshot='wisting_change_swat_27_layer_4.png')

# mesh = pv.PolyData(vertices)
# pl = pv.Plotter()
# pl.background_color = "white"
# pl.camera_position = 'xy'
# pl.camera.view_angle = 60.0
# orig_map = plt.cm.get_cmap('jet')
# reversed_map = orig_map.reversed()
# sargs = dict(height=0.25, vertical=True, position_x=0.85, position_y=0.05)
# pl.add_mesh(mesh, show_edges=True, cmap=reversed_map, scalars=cell_values, scalar_bar_args=sargs, reset_camera=True)
# pl.show()
