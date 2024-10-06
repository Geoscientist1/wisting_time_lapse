import numpy as np
import math as m
import scipy.integrate as integrate
from scipy.integrate import quad
from numpy import *
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))
import matplotlib.pyplot as plt
from matplotlib import ticker

# import cmocean

# f: means the finer model (314*318*101)
# c: means the coarser model (157*159*101)

model_resolution = "c"

if model_resolution == "c":
    A = 157
    B = 159
    C = 101
else:
    A = 314
    B = 318


def readfiles(f) -> object:
    with open(f, "r") as file:
        array = list()
        arr = []
        for line in file:
            array.append(line.split())
        for i in range(len(array)):
            for j in range(len(array[i])):
                c = array[i][j]
                s = []
                Cond = False
                for char in c:
                    if char == '*':
                        Cond = True
                        a, b = c.split('*')
                        a = int(a)
                        b = float(b)
                        b = int(b)
                        temp = [b for i in range(a)]
                        arr.append(temp)
                        break
                    else:
                        s.append(char)
                        f = ''.join(s)
                f = float(f)
                value = f
                array[i][j] = value
                for k in c:
                    if k != '*' and Cond == False:
                        arr.append(f)
                        break
        f_arr = np.array(arr, dtype=object)
        f_arr = f_arr.flatten()
        f_arr = np.hstack(f_arr)

    return f_arr


# data = np.loadtxt('SAT_27.txt', dtype='str')  # SAT_27.txt is the file containing the water saturation values
# SAT_27_f.txt: means the finer model (314*318*101)
# SAT_27_c.txt: means the coarser model (157*159*101)
f27 = 'SAT_27_' + model_resolution + '.txt'
f29 = 'SAT_29_' + model_resolution + '.txt'
f31 = 'SAT_31_' + model_resolution + '.txt'
f32 = 'SAT_32_' + model_resolution + '.txt'
f33 = 'SAT_33_' + model_resolution + '.txt'
f35 = 'SAT_35_' + model_resolution + '.txt'
f37 = 'SAT_37_' + model_resolution + '.txt'
f39 = 'SAT_39_' + model_resolution + '.txt'
f41 = 'SAT_41_' + model_resolution + '.txt'
f42 = 'SAT_42_' + model_resolution + '.txt'
f43 = 'SAT_43_' + model_resolution + '.txt'
f45 = 'SAT_45_' + model_resolution + '.txt'
f47 = 'SAT_47_' + model_resolution + '.txt'
f49 = 'SAT_49_' + model_resolution + '.txt'
f51 = 'SAT_51_' + model_resolution + '.txt'
f52 = 'SAT_52_' + model_resolution + '.txt'
f53 = 'SAT_53_' + model_resolution + '.txt'
f55 = 'SAT_55_' + model_resolution + '.txt'
f57 = 'SAT_57_' + model_resolution + '.txt'

# data = np.loadtxt('SATG_27.txt', dtype='str')  # SOIL_27.txt is the file containing the oil saturation values
so27 = 'SOIL_27_' + model_resolution + '.txt'
so29 = 'SOIL_29_' + model_resolution + '.txt'
so31 = 'SOIL_31_' + model_resolution + '.txt'
so32 = 'SOIL_32_' + model_resolution + '.txt'
so33 = 'SOIL_33_' + model_resolution + '.txt'
so35 = 'SOIL_35_' + model_resolution + '.txt'
so37 = 'SOIL_37_' + model_resolution + '.txt'
so39 = 'SOIL_39_' + model_resolution + '.txt'
so42 = 'SOIL_42_' + model_resolution + '.txt'
so45 = 'SOIL_45_' + model_resolution + '.txt'
so47 = 'SOIL_47_' + model_resolution + '.txt'
so49 = 'SOIL_49_' + model_resolution + '.txt'
so52 = 'SOIL_52_' + model_resolution + '.txt'
so53 = 'SOIL_53_' + model_resolution + '.txt'
so57 = 'SOIL_57_' + model_resolution + '.txt'

# data = np.loadtxt('SATG_27.txt', dtype='str')  # SATG_27.txt is the file containing the gas saturation values
fg27 = 'SGAS_27_' + model_resolution + '.txt'
fg29 = 'SGAS_29_' + model_resolution + '.txt'
fg31 = 'SGAS_31_' + model_resolution + '.txt'
fg32 = 'SGAS_32_' + model_resolution + '.txt'
fg33 = 'SGAS_33_' + model_resolution + '.txt'
fg35 = 'SGAS_35_' + model_resolution + '.txt'
fg37 = 'SGAS_37_' + model_resolution + '.txt'
fg39 = 'SGAS_39_' + model_resolution + '.txt'
fg42 = 'SGAS_42_' + model_resolution + '.txt'
fg45 = 'SGAS_45_' + model_resolution + '.txt'
fg47 = 'SGAS_47_' + model_resolution + '.txt'
fg49 = 'SGAS_49_' + model_resolution + '.txt'
fg52 = 'SGAS_52_' + model_resolution + '.txt'
fg53 = 'SGAS_53_' + model_resolution + '.txt'
fg57 = 'SGAS_57_' + model_resolution + '.txt'

vwater_26 = 'vwater_26.txt'  # Fluid in place water FIPWAT
vwater_27 = 'vwater_27.txt'
vwater_29 = 'vwater_29.txt'
vwater_31 = 'vwater_31.txt'
vwater_33 = 'vwater_33.txt'
vwater_35 = 'vwater_35.txt'
vwater_39 = 'vwater_39.txt'
vwater_45 = 'vwater_45.txt'
vwater_49 = 'vwater_49.txt'
vwater_53 = 'vwater_53.txt'
vwater_57 = 'vwater_57.txt'

voil_26 = 'VOIL_26.txt'  # Fluid in place oil FIPOIL
voil_27 = 'VOIL_27.txt'
voil_29 = 'VOIL_29.txt'
voil_31 = 'VOIL_31.txt'
voil_33 = 'VOIL_33.txt'
voil_35 = 'VOIL_35.txt'
voil_39 = 'VOIL_39.txt'
voil_45 = 'VOIL_45.txt'
voil_49 = 'VOIL_49.txt'
voil_53 = 'VOIL_53.txt'
voil_57 = 'VOIL_57.txt'
#
vgas_27 = 'vgas_27.txt'  # Fluid in place oil FIPGAS
vgas_29 = 'vgas_29.txt'
vgas_31 = 'vgas_31.txt'
vgas_33 = 'vgas_33.txt'
vgas_35 = 'vgas_35.txt'
vgas_39 = 'vgas_39.txt'
vgas_45 = 'vgas_45.txt'
vgas_49 = 'vgas_49.txt'
vgas_53 = 'vgas_53.txt'
vgas_57 = 'vgas_57.txt'

# data = np.loadtxt('DX.txt', dtype='str')  # DX.txt is the file containing the Cell X dimensions which are fixed
# throughout the years
f1 = 'DX_' + model_resolution + '.txt'
f4 = 'dx_accum_' + model_resolution + '.txt'

# data = np.loadtxt('DY.txt', dtype='str')  # DY.txt is the file containing the Cell Y dimensions which are fixed
# throughout the years
f2 = 'DY_' + model_resolution + '.txt'
f5 = 'dy_accum_' + model_resolution + '.txt'

# data = np.loadtxt('DZ.txt', dtype='str')  # DZ.txt is the file containing the Cell Z dimensions which are fixed
# throughout the years
f3 = 'DZ_' + model_resolution + '.txt'

# Pore volume at reference conditions
pv = 'PORV_' + model_resolution + '.txt'

# Cell top depth
celltopdpt = 'cell_top_dpt_' + model_resolution + '.txt'


# Cell center depth
cel_cen_dpt = 'cell_center_depth_' + model_resolution + '.txt'

# Depth_map of Stø(contours)
dpt = 'depth_map_' + model_resolution + '.txt'

# Cell dimensions dx, dy and dz
dx = readfiles(f1)
# np.savetxt('dx.txt', dx, delimiter=' ', header='dx data')
dy = readfiles(f2)
# np.savetxt('dy.txt', dy, delimiter=' ', header='dy data')
dz = readfiles(f3)
# np.savetxt('dz.txt', dz, delimiter=' ', header='dz data')
dx_accum = readfiles(f4)
dy_accum = readfiles(f5)
celltopdpt = readfiles(celltopdpt)
cel_cen_dpt = readfiles(cel_cen_dpt)
# Depth_map of Stø(contours)
dpt = readfiles(dpt)


# satXY contains water saturation values corresponding to year XY
sat27 = readfiles(f27)
# np.savetxt('sat27.txt', sat27, delimiter=' ', header='sat27 data')
sat29 = readfiles(f29)
# np.savetxt('sat29.txt', sat29, delimiter=' ', header='sat29 data')
sat31 = readfiles(f31)
# np.savetxt('sat31.txt', sat31, delimiter=' ', header='sat31 data')
sat32 = readfiles(f32)
# np.savetxt('sat32.txt', sat32, delimiter=' ', header='sat32 data')
sat33 = readfiles(f33)
# np.savetxt('sat33.txt', sat33, delimiter=' ', header='sat33 data')
sat35 = readfiles(f35)
# np.savetxt('sat35.txt', sat35, delimiter=' ', header='sat35 data')
sat37 = readfiles(f37)
# np.savetxt('sat37.txt', sat37, delimiter=' ', header='sat37 data')
sat39 = readfiles(f39)
# np.savetxt('sat39.txt', sat39, delimiter=' ', header='sat39 data')
sat41 = readfiles(f41)
# np.savetxt('sat41.txt', sat41, delimiter=' ', header='sat41 data')
sat42 = readfiles(f42)
# np.savetxt('sat42.txt', sat42, delimiter=' ', header='sat42 data')
sat43 = readfiles(f43)
# np.savetxt('sat43.txt', sat43, delimiter=' ', header='sat43 data')
sat45 = readfiles(f45)
# np.savetxt('sat45.txt', sat45, delimiter=' ', header='sat45 data')
sat47 = readfiles(f47)
# np.savetxt('sat47.txt', sat47, delimiter=' ', header='sat47 data')
sat49 = readfiles(f49)
# np.savetxt('sat49.txt', sat49, delimiter=' ', header='sat49 data')
sat51 = readfiles(f51)
# np.savetxt('sat51.txt', sat51, delimiter=' ', header='sat51 data')
sat52 = readfiles(f52)
# np.savetxt('sat52.txt', sat52, delimiter=' ', header='sat52 data')
sat53 = readfiles(f53)
# np.savetxt('sat53.txt', sat53, delimiter=' ', header='sat53 data')
sat55 = readfiles(f55)
# np.savetxt('sat55.txt', sat55, delimiter=' ', header='sat55 data')
sat57 = readfiles(f57)
# np.savetxt('sat57.txt', sat57, delimiter=' ', header='sat57 data')

# satXY contains oil saturation values corresponding to year XY
soil_27 = readfiles(so27)
# np.savetxt('soil_27.txt', soil_27, delimiter=' ', header='soil_27 data')
soil_29 = readfiles(so29)
# np.savetxt('soil_29.txt', soil_29, delimiter=' ', header='soil_29 data')
soil_31 = readfiles(so31)
# np.savetxt('soil_31.txt', soil_31, delimiter=' ', header='soil_31 data')
soil_32 = readfiles(so32)
# np.savetxt('soil_32.txt', soil_32, delimiter=' ', header='soil_32 data')
soil_33 = readfiles(so33)
# np.savetxt('soil_33.txt', soil_33, delimiter=' ', header='soil_33 data')
soil_35 = readfiles(so35)
# np.savetxt('soil_35.txt', soil_35, delimiter=' ', header='soil_35 data')
soil_37 = readfiles(so37)
# np.savetxt('soil_37.txt', soil_37, delimiter=' ', header='soil_37 data')
soil_39 = readfiles(so39)
# np.savetxt('soil_39.txt', soil_39, delimiter=' ', header='soil_39 data')
soil_42 = readfiles(so42)
# np.savetxt('soil_42.txt', soil_42, delimiter=' ', header='soil_42 data')
soil_45 = readfiles(so45)
# np.savetxt('soil_45.txt', soil_45, delimiter=' ', header='soil_45 data')
soil_47 = readfiles(so47)
# np.savetxt('soil_47.txt', soil_47, delimiter=' ', header='soil_47 data')
soil_49 = readfiles(so49)
# np.savetxt('soil_49.txt', soil_49, delimiter=' ', header='soil_49 data')
soil_52 = readfiles(so52)
# np.savetxt('soil_52.txt', soil_52, delimiter=' ', header='soil_52 data')
soil_53 = readfiles(so53)
# np.savetxt('soil_53.txt', soil_53, delimiter=' ', header='soil_53 data')
soil_57 = readfiles(so57)
# np.savetxt('soil_57.txt', soil_57, delimiter=' ', header='soil_57 data')

# satXY contains gas saturation values corresponding to year XY
sgas_27 = readfiles(fg27)
# np.savetxt('sgas_27.txt', sgas_27, delimiter=' ', header='sgas_27 data')
sgas_29 = readfiles(fg29)
# np.savetxt('sgas_29.txt', sgas_29, delimiter=' ', header='sgas_29 data')
sgas_31 = readfiles(fg31)
# np.savetxt('sgas_31.txt', sgas_31, delimiter=' ', header='sgas_31 data')
sgas_32 = readfiles(fg32)
# np.savetxt('sgas_32.txt', sgas_32, delimiter=' ', header='sgas_32 data')
sgas_33 = readfiles(fg33)
# np.savetxt('sgas_33.txt', sgas_33, delimiter=' ', header='sgas_33 data')
sgas_35 = readfiles(fg35)
# np.savetxt('sgas_35.txt', sgas_35, delimiter=' ', header='sgas_35 data')
sgas_37 = readfiles(fg37)
# np.savetxt('sgas_37.txt', sgas_37, delimiter=' ', header='sgas_37 data')
sgas_39 = readfiles(fg39)
# np.savetxt('sgas_39.txt', sgas_39, delimiter=' ', header='sgas_39 data')
sgas_42 = readfiles(fg42)
# np.savetxt('sgas_42.txt', sgas_42, delimiter=' ', header='sgas_42 data')
sgas_45 = readfiles(fg45)
# np.savetxt('sgas_45.txt', sgas_45, delimiter=' ', header='sgas_45 data')
sgas_47 = readfiles(fg47)
# np.savetxt('sgas_47.txt', sgas_47, delimiter=' ', header='sgas_47 data')
sgas_49 = readfiles(fg49)
# np.savetxt('sgas_49.txt', sgas_49, delimiter=' ', header='sgas_49 data')
sgas_52 = readfiles(fg52)
# np.savetxt('sgas_52.txt', sgas_52, delimiter=' ', header='sgas_52 data')
sgas_53 = readfiles(fg53)
# np.savetxt('sgas_53.txt', sgas_53, delimiter=' ', header='sgas_53 data')
sgas_57 = readfiles(fg57)
# np.savetxt('sgas_57.txt', sgas_57, delimiter=' ', header='sgas_57 data')

# # vlXY and vo_XY contain water and oil volume in place values respectively in each cell corresponding to year XY
vw_26 = readfiles(vwater_26)
vw_27 = readfiles(vwater_27)
vw_29 = readfiles(vwater_29)
vw_31 = readfiles(vwater_31)
vw_33 = readfiles(vwater_33)
vw_35 = readfiles(vwater_35)
vw_39 = readfiles(vwater_39)
vw_45 = readfiles(vwater_45)
vw_49 = readfiles(vwater_49)
vw_53 = readfiles(vwater_53)
vw_57 = readfiles(vwater_57)

vg_27 = readfiles(vgas_27)
vg_29 = readfiles(vgas_29)
vg_31 = readfiles(vgas_31)
vg_33 = readfiles(vgas_33)
vg_35 = readfiles(vgas_35)
vg_39 = readfiles(vgas_39)
vg_45 = readfiles(vgas_45)
vg_49 = readfiles(vgas_49)
vg_53 = readfiles(vgas_53)
vg_57 = readfiles(vgas_57)

vo_26 = readfiles(voil_26)
vo_27 = readfiles(voil_27)
vo_29 = readfiles(voil_29)
vo_31 = readfiles(voil_31)
vo_33 = readfiles(voil_33)
vo_35 = readfiles(voil_35)
vo_39 = readfiles(voil_39)
vo_45 = readfiles(voil_45)
vo_49 = readfiles(voil_49)
vo_53 = readfiles(voil_53)
vo_57 = readfiles(voil_57)
########################################################################################################################

# porv contains pore volume values in each cell
porv = readfiles(pv)

# Here the porosity in each cell is calculated
phi = np.zeros(len(porv))
cell_vol = dx * dy * dz  # Calculating the volume of each cell in order to find the porosity in each cell
for m in range(len(porv)):
    if cell_vol[m] == 0:
        phi[m] = 0
    else:
        phi[m] = porv[m] / cell_vol[m]

# np.savetxt('cellVol_' + model_resolution + '.txt', cell_vol)
# np.savetxt('porosity_' + model_resolution + '.txt', phi)


########################################################################################################################
# This function averages the input function
def avg(list1, nx, ny, nz):
    temp1 = 0
    temp2 = 0
    new = []
    new1 = []
    m = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                temp1 += list1[m]
                if np.mod(k, 9) == 0:
                    b = temp1 / 10
                    new.append(b)
                    temp1 = 0
                m += 1
    s = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(11):
                temp2 += new[s]
                if np.mod(k, 11) == 0:
                    new1.append(temp2 / 12)
                    temp2 = 0
                s += 1
    return np.array(new1)


def avg_poro(list1, nx, ny, nz):
    temp1 = 0
    temp2 = 0
    new = []
    new1 = []
    m = 0
    n = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if list1[m] != 0:
                    n += 1
                    temp1 += list1[m]
                if np.mod(n, 9) == 0:
                    b = temp1 / 10
                    new.append(b)
                    temp1 = 0
                m += 1
    count = []
    for i in range(len(new)):
        if new[i] == 0:
            count.append(i)
    np.array(count)
    new2 = np.delete(new, count)
    s = 0
    w = np.floor((len(new2)) / (157 * 159))
    if w == 1:
        for i in range(nx):
            for j in range(ny):
                for k in range(int(w)):
                    temp2 += new2[s]
                    new1.append(temp2 / w)
                    temp2 = 0
                    s += 1
    else:
        for i in range(nx):
            for j in range(ny):
                for k in range(int(w)):
                    temp2 += new2[s]
                    if np.mod(k, w - 1) == 0:
                        new1.append(temp2 / (w + 1))
                        temp2 = 0
                    s += 1
    return np.array(new1)


dx_avg = avg(dx, 157, 159, 101)
dy_avg = avg(dy, 157, 159, 101)
# np.savetxt('dy_avg.txt', dy_avg, delimiter=' ', header='dy_avg data')
dz_avg = avg(dz, 157, 159, 101)
# np.savetxt('dz_avg.txt', dz_avg, delimiter=' ', header='dz_avg data')
# poro_avg = avg_poro(poro, 157, 159, 101)
# np.savetxt('porosity_avg.txt', poro_avg)
sat27_avg = avg(sat27, 157, 159, 101)
sat29_avg = avg(sat29, 157, 159, 101)
sat31_avg = avg(sat31, 157, 159, 101)
sat33_avg = avg(sat33, 157, 159, 101)
sat35_avg = avg(sat35, 157, 159, 101)
sat39_avg = avg(sat39, 157, 159, 101)
sat45_avg = avg(sat45, 157, 159, 101)
sat49_avg = avg(sat49, 157, 159, 101)
sat53_avg = avg(sat53, 157, 159, 101)
sat57_avg = avg(sat57, 157, 159, 101)
########################################################################################################################
