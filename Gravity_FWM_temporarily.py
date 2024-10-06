
def mod_grid_coarsening(dt):
    dt_coarse = np.zeros((k_range_new, B, A))
    dt = np.reshape(dt, (k_range, B, A))
    for i in range(A):
        for j in range(B):
            s = 0
            fk = 0
            temp = 0
            for k in range(k_range):
                if dt[k, j, i] == 0:
                    summen = 0
                    count = 0
                    for n in range(k_range):
                        summen += dt[n, j, i]
                        if dt[n, j, i] != 0:
                            count += 1
                    if count != 0:
                        dt[k, j, i] = summen / count
                        X = dt[k, j, i]
                # if dt[k, j, i] != 0:
                #     fk += 1
                # temp += dt[k, j, i]
                # if k != 0 and mod(k, mod_fk) == 0 and fk != 0:
                #     Y = dt[k, j, i]
                #     dt_coarse[s, j, i] = temp / fk
                #     s += 1
                #     temp = 0
                #     fk = 0
    result = dt.flatten()
    return result





# dx_uni, dy_uni, dz_uni, dx_accum_uni, dy_accum_uni = uniform_horizontal_gridding()

# np.savetxt('dx_uniform.txt', dx_uni, delimiter=' ', header='Uniform DX values')
# np.savetxt('dy_uniform.txt', dy_uni, delimiter=' ', header='Uniform DY values')
# np.savetxt('dz_uniform.txt', dz_uni, delimiter=' ', header='Uniform DZ values')
# np.savetxt('dx_accum_uniform.txt', dx_accum_uni, delimiter=' ', header='Uniform DX_accum values')
# np.savetxt('dy_accum_uniform.txt', dy_accum_uni, delimiter=' ', header='Uniform Dy_accum values')
########################################################################################################################
# This is to choose if the grid is uniform (dx = dy = 70m, and dz = 2m)
# dx_coarse = mod_grid_coarsening(dx_uni)
# dy_coarse = mod_grid_coarsening(dy_uni)
# dz_coarse = mod_grid_coarsening(dz_uni)
# dx_accum_coarse = mod_grid_coarsening(dx_accum)
# dy_accum_coarse = mod_grid_coarsening(dy_accum)
########################################################################################################################
# This is to choose if the grid is non-uniform (Petrel grid values for dx, dy and dz)
dx_coarse = mod_grid_coarsening(dx)
dy_coarse = mod_grid_coarsening(dy)
dz_coarse = mod_grid_coarsening_dz(dz)
dx_accum_coarse = mod_grid_coarsening(dx_accum)
dy_accum_coarse = mod_grid_coarsening(dy_accum)
#######################################################################################################################
dv_coarse = mod_grid_coarsening(dv)
cel_cen_dpt_coarse = mod_grid_coarsening(cel_cen_dpt)
phi_coarse = mod_grid_coarsening(phi)
rho0_coarse = rho0_array_generator()
rho27_27_coarse = mod_grid_coarsening_rho(rho27_27)
rho27_29_coarse = mod_grid_coarsening_rho(rho27_29)
rho27_31_coarse = mod_grid_coarsening_rho(rho27_31)
rho27_33_coarse = mod_grid_coarsening_rho(rho27_33)
rho27_35_coarse = mod_grid_coarsening_rho(rho27_35)
rho27_39_coarse = mod_grid_coarsening_rho(rho27_39)
rho27_45_coarse = mod_grid_coarsening_rho(rho27_45)
rho27_49_coarse = mod_grid_coarsening_rho(rho27_49)
rho27_53_coarse = mod_grid_coarsening_rho(rho27_53)
rho27_57_coarse = mod_grid_coarsening_rho(rho27_57)

np.savetxt('dx_coarse.txt', dx_coarse, delimiter=' ', header='')
np.savetxt('dy_coarse.txt', dy_coarse, delimiter=' ', header='')
np.savetxt('dz_coarse.txt', dz_coarse, delimiter=' ', header='')
np.savetxt('dv_coarse.txt', dv_coarse, delimiter=' ', header='')
np.savetxt('dx_accum_coarse.txt', dx_accum_coarse, delimiter=' ', header='')
np.savetxt('dy_accum_coarse.txt', dy_accum_coarse, delimiter=' ', header='')
np.savetxt('cel_cen_dpt_coarse.txt', cel_cen_dpt_coarse, delimiter=' ', header='')
np.savetxt('phi_coarse.txt', phi_coarse, delimiter=' ', header='')


# Dot not forget to generate a density array consisting of rho0 only, the code for that is now written
# Do not forget to import the density arrays to Grav_Modeling_coarse_grid script and use the in the FWM engine.
# Do not forget to test the code developed here for FWM for the same densisties without rho0 just to see that the same
# results are reproducible.
# Then generate the delta_g for the different densities and for the rho0
# Then substract the two maps from each other, this should be the result.


# Gravity_Modeling_coarse_grid modified script from the original one
















sat27 = np.loadtxt('sat27_coarse.txt', delimiter=' ')
sat29 = np.loadtxt('sat29_coarse.txt', delimiter=' ')
sat31 = np.loadtxt('sat31_coarse.txt', delimiter=' ')
sat33 = np.loadtxt('sat33_coarse.txt', delimiter=' ')
sat35 = np.loadtxt('sat35_coarse.txt', delimiter=' ')
sat39 = np.loadtxt('sat39_coarse.txt', delimiter=' ')
sat45 = np.loadtxt('sat45_coarse.txt', delimiter=' ')
sat49 = np.loadtxt('sat49_coarse.txt', delimiter=' ')
sat53 = np.loadtxt('sat53_coarse.txt', delimiter=' ')
sat57 = np.loadtxt('sat57_coarse.txt', delimiter=' ')

soil_27 = np.loadtxt('soil_27_coarse.txt', delimiter=' ')
soil_29 = np.loadtxt('soil_29_coarse.txt', delimiter=' ')
soil_31 = np.loadtxt('soil_31_coarse.txt', delimiter=' ')
soil_33 = np.loadtxt('soil_33_coarse.txt', delimiter=' ')
soil_35 = np.loadtxt('soil_35_coarse.txt', delimiter=' ')
soil_39 = np.loadtxt('soil_39_coarse.txt', delimiter=' ')
soil_45 = np.loadtxt('soil_45_coarse.txt', delimiter=' ')
soil_49 = np.loadtxt('soil_49_coarse.txt', delimiter=' ')
soil_53 = np.loadtxt('soil_53_coarse.txt', delimiter=' ')
soil_57 = np.loadtxt('soil_57_coarse.txt', delimiter=' ')

rho_oil = 835  # Oil density
rho_water = 1040  # Water density
d_rho = rho_water - rho_oil  # Density difference between water and oil
depth = 230  # Reservoir depth at Wisting
v = 70 * 1e6  # Total volume of reserves corresponding to 440 million barrels of oil
G = 6.67 * 1e-11  # Gravity constant
k_range = 101
k_range_new = 10


# Now we account for the different water saturations corresponding to the different times from year 2027 all the way
# to year 2039
################################################################################
################################################################################
def grav_mod_coarse(rho):
    XSCALE = 84
    YSCALE = 89

    # sat_1 = np.reshape(sat_1, (k_range_new, B, A))
    # sat_2 = np.reshape(sat_2, (k_range_new, B, A))
    # soil_1 = np.reshape(soil_1, (k_range_new, B, A))
    # soil_2 = np.reshape(soil_2, (k_range_new, B, A))
    # cel_cen_dpt_mo = np.reshape(cel_cen_dpt, (k_range_new, B, A))
    # phi_mo = np.reshape(phi, (k_range_new, B, A))
    # dx_mo = np.reshape(dx, (k_range_new, B, A))
    # dy_mo = np.reshape(dy, (k_range_new, B, A))
    # dz_mo = np.reshape(dz, (k_range_new, B, A))
    # dx_accum_mo = np.reshape(dx_accum, (k_range_new, B, A))
    # dy_accum_mo = np.reshape(dy_accum, (k_range_new, B, A))

    delta_g = np.zeros((B, A))
    C = 0
    rho_1 = ((1 - soil_1) * rho_water) + ((1 - sat_1) * rho_oil)
    rho_2 = ((1 - soil_2) * rho_water) + ((1 - sat_2) * rho_oil)
    delta_rho = rho_2 - rho_1
    a = G * (cel_cen_dpt - 400) * rho * phi
    b = dx * dy * dz  # this corresponds to V v(volume) in the formula for calculating delta gz
    c = a * b
    for y_step in range(0, B):
        for x_step in range(0, A):
            print("y", y_step)
            print("x", x_step)
            m = 0
            for k in range(k_range):
                for j in range(B):
                    for i in range(A):
                        r = np.sqrt(((x_step * XSCALE) - dx_accum[m]) ** 2 + ((y_step * YSCALE) - dy_accum[m]) ** 2)
                        d = pow((pow(r, 2) + pow(cel_cen_dpt[m] - 400, 2)), 1.5)
                        delta_g[158 - y_step, x_step] += ((c[m]) / d)
                        C += ((c[m]) / d)
                        m += 1
    return delta_g * 1e8


# gv27_27 = grav_mod_coarse(sat27, sat27, soil_27, soil_27)
# gv27_29 = grav_mod_coarse(sat27, sat29, soil_27, soil_29)
# gv27_31 = grav_mod_coarse(sat27, sat31, soil_27, soil_31)
# gv27_33 = grav_mod_coarse(sat27, sat33, soil_27, soil_33)rigin
# gv27_35 = grav_mod_coarse(sat27, sat35, soil_27, soil_35)
# gv27_39 = grav_mod_coarse(sat27, sat39, soil_27, soil_39)
# gv27_45 = grav_mod_coarse(sat27, sat45, soil_27, soil_45)
# gv27_49 = grav_mod_coarse(sat27, sat49, soil_27, soil_49)
# gv27_53 = grav_mod_coarse(sat27, sat53, soil_27, soil_53)
gv27_57 = grav_mod_coarse(rho27_57)
gv_rho0 = grav_mod_coarse(rho0)
gv_main_result = gv27_57 - gv_rho0

np.savetxt('grav_data_57_coarse_non_uniform_real_101_avg_2_cells_full_model.txt', gv_main_result, delimiter=' ',
           header='Gravity data Coarse Grid for year 2057')
