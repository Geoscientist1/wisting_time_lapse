# Load data

import plotly.graph_objects as go
import dash
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

df = pd.read_csv('myOutFile_gslib2sblwiz_overburden.txt', sep=' ')
app = dash.Dash(__name__)

x = df.XA
y = df.YA
z = df.ZA
g = df.GA

# Simple Scatter
plt.scatter(x, z)
plt.show()

# create meshgrid
xi = np.arange(x.min(), x.max(), (x.max() - x.min()) / 100)
yi = np.arange(y.min(), y.max(), (y.max() - y.min()) / 100)
xi, yi = np.meshgrid(xi, yi)
# interpolate
zi = griddata((x, y), z, (xi, yi), method='linear')
########################################################################################################################





########################################################################################################################

# Plot data
fig = go.Figure(data=[
    go.Surface(z=zi, x=xi, y=yi, opacity=0.75, colorscale='Jet'),
])

fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))

zoom = 1.25
fig.update_layout(hovermode='closest', xaxis=dict(showspikes=False),
                  scene_camera=dict(eye=dict(x=0.5 * zoom, y=-1 * zoom, z=0.5 * zoom)))
fig.show()

# if __name__ == '__main__':
#     app.run_server(debug=True, use_reloader=False, port=8051)
