import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as plc
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.ndimage import gaussian_filter

# create some fake data
x = np.linspace(0, 13213.48, 157)
y = np.linspace(0, 14173.40, 159)
sat27 = np.loadtxt('sat27_coarse.txt', delimiter=' ')
sat57 = np.loadtxt('sat57_coarse.txt', delimiter=' ')
ref_model = np.loadtxt('reference_model.txt', delimiter=' ')
# dt = np.loadtxt('recovered_model_temp_118_9_10000.txt', delimiter=' ')
dt = np.loadtxt("start_model_146.txt", delimiter=" ")
dswat_observed = np.subtract(sat57, sat27)
dswat_recovered = np.subtract(dt, sat27)
# dswat_recovered = dt
# dt = dswat_observed
# here are the x,y and respective z values
X, Y = np.meshgrid(x, y)
# zrng = np.loadtxt('cel_cen_dpt_coarse.txt', delimiter=' ') # This is used for the real model
zrng = np.loadtxt('dpt_cell_initial_model_inversion.txt', delimiter=' ')  # This is used for the start model
Z = np.zeros((5, 159, 157))
m = 0
for k in range(5):
    for j in range(159):
        for i in range(157):
            if zrng[m] == 1:
                Z[k, j, i] = None
                m += 1
            else:
                Z[k, j, i] = -zrng[m]
                m += 1
# this is the value to use for the color
V = np.reshape(dt, (5, 159, 157))
# create the figure, add a 3d axis, set the viewing angle
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(70, 60)

# here we create the surface plot, but pass V through a colormap
# to create a different color for each patch

maxlimit = 0.5
minlimit = 0.0

for i in range(5):
    norm = plc.Normalize(vmin=minlimit, vmax=maxlimit)
    surface = ax.plot_surface(X * 0.001, Y * 0.001, Z[i], facecolors=cm.jet(norm(V[i])), cmap=cm.jet, linewidth=0,
                              antialiased=False, vmin=minlimit, vmax=maxlimit)

# layer = 2
# norm = plc.Normalize(vmin=minlimit, vmax=maxlimit)
# surface = ax.plot_surface(X*0.001, Y*0.001, Z[layer], facecolors=cm.jet(norm(V[layer])), cmap=cm.jet, linewidth=0,
#                               antialiased=False, vmin=minlimit, vmax=maxlimit)

ax.set(zlim=(-1000, -500), zlabel="Depth [m]")
surface.set_clim(minlimit, maxlimit)
m = cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
fig.colorbar(m, shrink=0.5, aspect=5)
plt.xlabel("East [km]")
ax.invert_yaxis()
plt.ylabel("North [km]")
ax.grid(False)
plt.show()
