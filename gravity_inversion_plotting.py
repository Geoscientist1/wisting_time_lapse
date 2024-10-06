import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import scipy.interpolate as si
import pandas as pd
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
import plotly.graph_objects as go
import plotly.express as px
import dash
import time
# from Grav_Modeling_coarse_grid import *
from gravity_inversion import *
from plotly.offline import plot  # for IDE use
app = dash.Dash(__name__)


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


########################################################################################################################
########################################################################################################################
def multiple3dsurfaceplots(model):
    model_reshaped = np.reshape(model, (z_range, y_range, x_range))
    X = np.linspace(0, x_max, x_range)
    Y = np.linspace(0, y_max, y_range)
    # Z = np.linspace(z_max, -300, z_range)
    Z = np.linspace(-300, z_max, z_range)
    model_to_be_drawn = np.zeros((z_range, y_range, x_range))
    x = np.zeros((z_range, y_range, x_range))
    y = np.zeros((z_range, y_range, x_range))
    z = np.zeros((z_range, y_range, x_range))

    # 20 cells for the water layer with uniform resistivity being equal to 0.1 ohm-m
    for k in range(0, z_range):
        for j in range(y_range):
            for i in range(x_range):
                model_to_be_drawn[k, j, i] = model_reshaped[k, j, i]
                z[k, j, i] = Z[k]
                x[k, j, i] = X[i]
                y[k, j, i] = Y[j]

    return model_to_be_drawn, x, y, z


########################################################################################################################
########################################################################################################################
# rho_starting_model = np.loadtxt("rho_initial.txt", delimiter=" ")
# rho_observed_model = np.loadtxt("rho_observed.txt", delimiter=" ")
# dt_initial = np.loadtxt("grav_inv_dummy_model_initial.txt", delimiter=" ")
# dt_observed = np.loadtxt("grav_inv_dummy_model_observed.txt", delimiter=" ")
########################################################################################################################
########################################################################################################################
def gravity_results_plotting_execution(rho_starting_model, rho_observed, recovered_model, recovered_data, data_observed):
    recovered_model, x_coordinates, y_coordinates, z_coordinates = multiple3dsurfaceplots(recovered_model)
    observed_model, x_coordinates, y_coordinates, z_coordinates = multiple3dsurfaceplots(rho_observed)
    initial_model, x_coordinates, y_coordinates, z_coordinates = multiple3dsurfaceplots(rho_starting_model)
    X = np.linspace(0, x_max, x_range)
    Y = np.linspace(0, y_max, y_range)
    Z = np.linspace(-1200, -300, 3)
    ####################################################################################################################
    ####################################################################################################################
    X = np.linspace(0, x_max, x_range)
    Y = np.linspace(0, y_max, y_range)
    plt.style.use('classic')

    ####################################################################################################################
    # diff = np.subtract(data_0, data_observed)
    # diff_noise = np.subtract(noise_data_observed, data_observed)
    # dt_labels_tmstp = [data_0, dt_observed, diff, noise_data_observed, diff_noise]
    # labels_tmstp = ["initial data", "observed data", "difference", "observed data with noise added",
    #                 "difference when noise included"]
    ####################################################################################################################
    diff = np.subtract(recovered_data, data_observed)
    dt_labels_tmstp = [data_observed, recovered_data, diff]
    labels_tmstp = ["observed data", "recovered_data", "difference"]
    labels_tmstp_0 = ["observed model", "recovered model"]
    ####################################################################################################################
    for i in range(3):
        dt = dt_labels_tmstp[i]
        fig = plt.figure(figsize=(6.4, 4.8))
        ax = fig.add_subplot(111)
        ax.set_aspect('auto')
        # f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
        # Z = f(Y * 0.001, X * 0.001)
        Z = gaussian_filter(dt, sigma=0, mode='reflect')
        plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
                cmap=plt.cm.get_cmap('jet', 1000),
                # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
                interpolation='nearest', origin='upper')
        clb = plt.colorbar()
        tick_locator = ticker.MaxNLocator(nbins=7)
        clb.locator = tick_locator
        clb.update_ticks()
        clb.ax.set_title('$(\u03BC$Gal)', fontsize=15)
        ax.set_ylabel("Y(m)", labelpad=15)
        if labels_tmstp[i] == "recovered_datad":
            print("")
            # plt.clim(0.04, 0.08)
        elif labels_tmstp[i] == "difference":
            print("")
            # plt.clim(0.04, 0.08)
        else:
            plt.clim(0.05, 0.27)
        # plt.axis('off')
        # plt.title(labels_tmstp[i] + "_layerstack_" + phs)
        plt.title(str(labels_tmstp[i]))
        plt.xlabel("X [$km$]")
        # ax.xaxis.tick_top()
        plt.ylabel("Y [$km$]")
        ax = plt.gca()
        ax.patch.set_visible(False)
        fig.patch.set_visible(False)
        # legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
        # ax.legend()
        # Put a nicer background color on the legend.
        # legend.get_frame()
        # ax.set_xticks(X)
        # ax.set_yticks(Y)
        # ax.grid(color='b', linestyle='-', linewidth=0.5)
        # ax.set_xticklabels(np.arange(1, 158, 1))
        # ax.set_yticklabels(np.arange(1, 160, 1))
        plt.show()
        fig.savefig(labels_tmstp[i], transparent=True)
    ####################################################################################################################
    ####################################################################################################################
    #   This is to be used if only the volume plot without slices should be displayed ##################################
    ####################################################################################################################
    name = 'eye = (x:2, y:2, z:0.1)'
    camera = dict(
        eye=dict(x=2, y=2, z=0.1)
    )
    density_models = [observed_model, recovered_model]
    plot_title_labels = ["3D Real Density Distribution", "3D Recovered Density Distribution"]
    for i in range(2):
        time.sleep(1)
        dt = density_models[i]
        dt_Smooth = gaussian_filter(dt, sigma=0, mode='reflect')
        fig = go.Figure(data=[
            go.Volume(x=x_coordinates.flatten(), y=y_coordinates.flatten(), z=z_coordinates.flatten(),
                value=dt_Smooth.flatten(), cmin=0, cmax=50, colorscale='Jet', opacity=.5, colorbar_len=0.3),
        ])
        fig.update_layout(scene_camera=camera, title=plot_title_labels[i], width=700, height=700)
        plot(fig, auto_open=True)
        fig.write_image(labels_tmstp_0[i] + ".png")
        layout = go.Layout(
        )
########################################################################################################################
########################################################################################################################
#   This is to be used if slices in volumetric data are wanted #########################################################
########################################################################################################################
# l = len(plotted_model) - 1
# plotted_model_Smooth = gaussian_filter(plotted_model, sigma=0.7, mode='reflect')
# fig = go.Figure(data=[
#     go.Volume(x=x_coordinates.flatten(), y=y_coordinates.flatten(), z=z_coordinates.flatten(),
#               value=plotted_model_Smooth.flatten(), cmin=0, cmax=140, colorscale='Jet', opacity=1, colorbar_len=0.8),
# ])
# frames = [go.Frame(data=go.Volume(
#     z=z_coordinates[l - k].flatten(), value=plotted_model_Smooth[l - k].flatten()),
#     name=str(k)) for k in range(len(plotted_model))]
# updatemenus = [dict(
#     buttons=[
#         dict(
#             args=[None, {"frame": {"duration": 100, "redraw": True},
#                          "fromcurrent": True, "transition": {"duration": 0}}],
#             label="Play",
#             method="animate"
#         ),
#         dict(
#             args=[[None], {"frame": {"duration": 0, "redraw": False},
#                            "mode": "immediate",
#                            "transition": {"duration": 0}}],
#             label="Pause",
#             method="animate"
#         )
#     ],
#     direction="left",
#     pad={"r": 10, "t": 87},
#     showactive=False,
#     type="buttons",
#     x=0.21,
#     xanchor="right",
#     y=-0.075,
#     yanchor="top"
# )]
#
# sliders = [dict(steps=[dict(method='animate',
#                             args=[[f'{k}'],
#                                   dict(mode='immediate',
#                                        frame=dict(duration=100, redraw=True),
#                                        transition=dict(duration=0))
#                                   ],
#                             label=f'{k + 1}'
#                             ) for k in range(len(plotted_model))],
#                 active=0,
#                 transition=dict(duration=0),
#                 x=0,  # slider starting position
#                 y=0,
#                 currentvalue=dict(font=dict(size=12),
#                                   prefix='frame: ',
#                                   visible=True,
#                                   xanchor='center'
#                                   ),
#                 len=1.0)  # slider length
#            ]
#
# fig.update_layout(title="3D Density Distribution",width=700, height=700, updatemenus=updatemenus, sliders=sliders)
# fig.update_layout(title="3D Density Distribution", width=700, height=700)
# fig.update(frames=frames)
# plot(fig, auto_open=True)
#
# layout = go.Layout(
#
# )
########################################################################################################################
########################################################################################################################
########################################################################################################################
# Define frames
# import plotly.graph_objects as go
#
# nb_frames = 2
#
# fig = go.Figure(frames=[go.Frame(data=go.Surface(
#     z=z_coordinates[k, :, :],
#     surfacecolor=plotted_model[k, :, :],
#     cmin=0, cmax=140
#     ),
#     name=str(k) # you need to name the frame for the animation to behave properly
#     )
#     for k in range(len(plotted_model))])
#
# # Add data to be displayed before animation starts
# fig.add_trace(go.Surface(
#     z=z_coordinates[2, :, :],
#     surfacecolor=plotted_model[2, :, :],
#     colorscale='Jet',
#     cmin=0, cmax=140,
#     colorbar=dict(thickness=20, ticklen=4)
#     ))
#
#
# def frame_args(duration):
#     return {
#             "frame": {"duration": duration},
#             "mode": "immediate",
#             "fromcurrent": True,
#             "transition": {"duration": duration, "easing": "linear"},
#         }
#
# sliders = [
#             {
#                 "pad": {"b": 10, "t": 60},
#                 "len": 0.9,
#                 "x": 0.1,
#                 "y": 0,
#                 "steps": [
#                     {
#                         "args": [[f.name], frame_args(0)],
#                         "label": str(k),
#                         "method": "animate",
#                     }
#                     for k, f in enumerate(fig.frames)
#                 ],
#             }
#         ]
#
# # Layout
# fig.update_layout(
#          title='Slices in volumetric data',
#          width=600,
#          height=600,
#          scene=dict(
#                     zaxis=dict(range=[-0.1, 6.8], autorange=False),
#                     aspectratio=dict(x=1, y=1, z=1),
#                     ),
#          updatemenus = [
#             {
#                 "buttons": [
#                     {
#                         "args": [None, frame_args(50)],
#                         "label": "&#9654;", # play symbol
#                         "method": "animate",
#                     },
#                     {
#                         "args": [[None], frame_args(0)],
#                         "label": "&#9724;", # pause symbol
#                         "method": "animate",
#                     },
#                 ],
#                 "direction": "left",
#                 "pad": {"r": 10, "t": 70},
#                 "type": "buttons",
#                 "x": 0.1,
#                 "y": 0,
#             }
#          ],
#          sliders=sliders
# )
#
# fig.show()
