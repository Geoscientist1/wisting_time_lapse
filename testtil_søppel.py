import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from tempfile import TemporaryFile
from scipy.ndimage import gaussian_filter
from matplotlib import ticker
from scipy.interpolate import RectBivariateSpline
from dsw_dpt_slices import GRID, POR, dz, k_range, fluid_type, typo
from files_readme import celltopdpt, phi, cell_vol
from files_readme import sat27, sat29, sat31, sat33, sat35, sat39, sat45, sat49, sat53, sat57
from files_readme import soil_27, soil_29, soil_31, soil_33, soil_35, soil_39, soil_45, soil_49, soil_53, soil_57
from files_readme import sgas_27, sgas_29, sgas_31, sgas_33, sgas_35, sgas_39, sgas_45, sgas_49, sgas_53, sgas_57
from files_readme import phi, cel_cen_dpt, dx, dy, dz, dx_accum, dy_accum

rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant


# Now we account for the different water saturations corresponding to the different times from year 2027 all the way
# to year 2039
################################################################################
################################################################################

# The function calculating the gravity responses delta gz
def grav_mod(sat_1, sat_2, soil_1, soil_2):
    M = 157
    N = 159
    XSCALE = 84
    YSCALE = 89

    delta_g = np.zeros((N, M))
    fixed1, fixed2, fixed3 = None, None, None
    rho_1 = ((1 - soil_1) * rho_water) + ((1 - sat_1) * rho_oil)
    rho_2 = ((1 - soil_2) * rho_water) + ((1 - sat_2) * rho_oil)
    delta_rho = rho_2 - rho_1
    a = G * (cel_cen_dpt - 400) * delta_rho * phi
    b = dx * dy * dz  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    for y_step in range(0, N):
        for x_step in range(0, M):
            print("y", y_step)
            print("x", x_step)
            m = 0
            for k in range(
                    1):  #########################################################################################
                for j in range(159):
                    for i in range(157):
                        r = np.sqrt(((x_step * XSCALE) - dx_accum[m]) ** 2 + ((y_step * YSCALE) - dy_accum[m]) ** 2)
                        d = pow((pow(r, 2) + pow(cel_cen_dpt[m] - 400, 2)), 1.5)
                        delta_g[157 - y_step, x_step] += ((c[m]) / d)
                        m += 1

            if x_step * XSCALE == 4284 and y_step * YSCALE == 12193:
                fixed1 = delta_g[157 - y_step, x_step]
            if x_step * XSCALE == 7980 and y_step * YSCALE == 6319:
                fixed2 = delta_g[157 - y_step, x_step]
            if x_step * XSCALE == 9996 and y_step * YSCALE == 8277:
                fixed3 = delta_g[157 - y_step, x_step]
    return delta_g * 1e8, fixed1 * 1e8, fixed2 * 1e8, fixed3 * 1e8


################################################################################
gv27_27, fixed_1_27, fixed_2_27, fixed_3_27 = grav_mod(sat27, sat27, soil_27, soil_27)
gv27_57, fixed_1_57, fixed_2_57, fixed_3_57 = grav_mod(sat27, sat57, soil_27, soil_57)

fixed_27 = np.array([fixed_1_27, fixed_2_27, fixed_3_27])
fixed_57 = np.array([fixed_1_57, fixed_2_57, fixed_3_57])


np.savetxt('søp_data_27.txt', gv27_27, delimiter=' ', header='Gravity data for year 2027')
np.savetxt('søp_data_57.txt', gv27_57, delimiter=' ', header='Gravity data for year 2057')

np.savetxt('søpfixed_27.txt', fixed_27, delimiter=' ', header='Station data for year 2027')

np.savetxt('søpfixed_57.txt', fixed_57, delimiter=' ', header='Station data for year 2057')


diff27_27 = gv27_27 - gv27_27
diff27_57 = gv27_57 - gv27_27
################################################################################
# gv27_31 = grav_mod(sat27, sat31)
# gv27_33 = grav_mod(sat27, sat33)
# gv27_35 = grav_mod(sat27, sat35)
# gv27_45 = grav_mod(sat27, sat45)
# gv27_49 = grav_mod(sat27, sat49)
# gv27_53 = grav_mod(sat27, sat53)
# print("min27-57", np.amin(gv27_57))
# print("max27-57", np.amax(gv27_57))
# gv27 = np.reshape(gv27, (157, 159))
# print(gv27_57)################################################################################
# Plotting of the gravity modeling results
# fig = plt.figure()
# dx = linspace(-50, 50, num=100, endpoint=True)
# dy = linspace(-50, 50, num=100, endpoint=True)
# plt.plot(dx, dy)

# sns.set(style="white")
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111)
# nptsx, nptsy = 100, 100
# dxg, dyg = np.meshgrid(np.linspace(dx.min(), dx.max(), nptsx), np.linspace(dy.min(), dy.max(), nptsy))
# triangles = tri.Triangulation(dx, dy)
# tri_interp = tri.CubicTriInterpolator(triangles, gv27)
# gv27g = tri_interp(dxg, dyg)
# change levels here according to your data# fig = plt.figure()
# levels = np.linspace(0, 10, 5)
# colormap = ax.contourf(dxg, dyg, gv27, levels,
#                      cmap=plt.cm.Blues,
#                     norm=plt.Normalize(vmax=gv27.max(), vmin=gv27.min()))
# plot data points
# ax.plot(dx, dy, color="#444444", marker="o", linestyle="", markersize=10)

# add a colorbar
# fig.colorbar(colormap, orientation='vertical', shrink=0.85)  # horizontal colour bar

# graph extras: look at xlim and ylim
# ax.set_xlim((0, 10))
# ax.set_ylim((0, 10))
# ax.set_aspect("equal", "box")
# plt.show()

# plt.ylabel('$y$', fontsize=18)
# #fig.savefig('')
# plt.colorbar()
# color_map = plt.imshow(gv27)

plt.style.use('classic')
X = np.linspace(0, 13213.48, 157)
Y = np.linspace(0, 14173.40, 159)
XSCALE = 84
YSCALE = 89

dt_labels_tmstp = [gv27_27, gv27_57]

diff_dt_labels_tmstp = [diff27_27, diff27_57]

labels_tmstp = ["gv_2727", "gv_2757"]

for i in range(2):
    dt = dt_labels_tmstp[i]
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')
    f = RectBivariateSpline(Y * 0.001, X * 0.001, dt)
    plt.plot(4284 * 0.001, 12193 * 0.001, marker='s', markersize=7, color="w", linestyle='None', label="Station 1")
    plt.plot(7980 * 0.001, 6319 * 0.001, marker='o', markersize=7, color="w", linestyle='None', label="Station 2")
    plt.plot(9996 * 0.001, 8277 * 0.001, marker='^', markersize=7, color="w", linestyle='None', label="Station 3")
    Z = f(Y * 0.001, X * 0.001)
    Z = gaussian_filter(dt, sigma=0, mode='reflect')
    plt.imshow(Z, extent=[min(X * 0.001), max(X * 0.001), min(Y * 0.001), max(Y * 0.001)],
               cmap=plt.cm.get_cmap('jet', 1000),
               # norm=LogNorm(vmin=max(dt[j].min(), LOGMIN)),
               interpolation='nearest', origin='upper')

    clb = plt.colorbar()
    tick_locator = ticker.MaxNLocator(nbins=7)
    clb.locator = tick_locator
    clb.update_ticks()
    # clb.set_label(r'log10(ATR)[$\Omega m^{2}$]', rotation=90, fontsize=18)
    clb.ax.set_title('$(\u03BC$Gal)', fontsize=15)
    ax.set_ylabel("Y(m)", labelpad=15)
    plt.clim(0, 1)  #####################################################################################
    # plt.title(labels_tmstp[i] + "_layerstack_" + phs)
    plt.xlabel("X [$km$]")
    # ax.xaxis.tick_top()
    plt.ylabel("Y [$km$]")
    ax = plt.gca()
    legend = ax.legend(loc='upper right', shadow=False, fontsize='x-small', numpoints=1)
    # ax.legend()
    # Put a nicer background color on the legend.
    # legend.get_frame()
    # ax.set_xticks(X)
    # ax.set_yticks(Y)
    ax.grid(color='b', linestyle='-', linewidth=0.5)
    # ax.set_xticklabels(np.arange(1, 158, 1))
    # ax.set_yticklabels(np.arange(1, 160, 1))
    plt.title(labels_tmstp[i])
    plt.show()
    fig.savefig(labels_tmstp[i])

    # legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    # plt.legend(loc="upper right")
    # plt.legend(loc=right, prop={'size': 6})
    # Put a nicer background color on the legend.
    # legend.get_frame()

diff = np.zeros(2)
diff[0] = sum(gv27_27 - gv27_27)
diff[1] = sum(gv27_57 - gv27_27)
print(diff)

time = [2027, 2057]
time1 = [27, 57]

# Bar chart showing the gravity changes from one timestep to the next
fig10 = plt.figure(figsize=(6.4, 4.8))
plt.grid()
x = np.arange(1, len(diff) + 1)
plt.bar(x, diff, align='center')
plt.xticks(x, time)
plt.xlabel('$Timesteps (years)$', fontsize=13)
plt.ylabel('$delta\_g\ (\u03BCGal)$', fontsize=13)
plt.show()
fig10.savefig("overview.pdf")

result1 = [fixed_1_27, fixed_1_57]
result2 = [fixed_2_27, fixed_2_57]
result3 = [fixed_3_27, fixed_3_57]

fig11, ax = plt.subplots()
x = np.arange(5)
ax.plot(result1, 'k--', label='Station 1')
ax.plot(result2, 'k:', label='Station 2')
ax.plot(result3, 'k', label='Station 3')
legend = ax.legend(loc='upper right', shadow=True, fontsize='small')
# Put a nicer background color on the legend.
legend.get_frame()
# plt.title('Gravity changes for three measuring stations')
plt.xlabel('$Timesteps (years)$', fontsize=18)
plt.ylabel('$delta\_g\ (\u03BCGal)$', fontsize=18)
plt.xticks(x, time)
plt.show()
fig11.savefig("final.pdf")

print('Result1', result1)
print('Result2', result2)
print('Result3', result3)
