import matplotlib.pyplot as plt
import numpy as np
from numpy import *
# import cv2
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
from Grav_Modeling_coarse_grid import grav_mod_coarse
from gravity_inversion import *
import time

start_time = time.time()
########################################################################################################################
########################################################################################################################
#         JACOBIAN CALCULATIONS
########################################################################################################################
########################################################################################################################
#         THE FIRST METHOD IS THE MODEL PERTURBATION METHOD
########################################################################################################################
print("WELCOME TO THE JACOBIAN CALCULATION SCRIPT, FROM HERE YOU WILL BE PROVIDED WITH THE JACOBIANS FOUND BY"
      "USING TWO DIFFERENT METHODS")


########################################################################################################################
########################################################################################################################
def jacobian_md_pert_method():
    N = x_range * y_range
    M = x_range * y_range * z_range
    jacobian = np.zeros((N, M))
    data_final_reshaped = data_final.reshape([N])

    for j in range(N):
        for i in range(M):
            print("This is the iteration number in terms of i:" + str(i) + " and in terms of j:" + str(j))
            perturbation = rho_final[i] * 0.01  # 1% perturbation when calculating the Jacobian
            temp = rho_final[i] + perturbation
            dt = grav_mod_coarse(temp, cell_dpt_starting_model, phi_starting_model, dx_starting_model,
                                 dy_starting_model, dz_starting_model, dx_accum_starting_model, dy_accum_starting_model,
                                 x_range, y_range, z_range, x_scale_starting_model, y_scale_starting_model)
            # The line just above is the FWM engine
            dt_reshaped = dt.reshape([N])
            approx = (data_final_reshaped[j] - dt_reshaped[j]) / (rho_final[i] - temp)
            jacobian[j, i] = approx
    return jacobian


########################################################################################################################
########################################################################################################################
jacobian_1 = jacobian_md_pert_method()
np.savetxt('jacobian_perturbation_method.txt', jacobian_1, delimiter=' ', header=' ')
########################################################################################################################
########################################################################################################################
end_time = time.time()
execution_time = end_time - start_time
print("EXECUTION TIME FOR THE JACOBIAN:", execution_time)
