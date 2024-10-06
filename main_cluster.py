from gravity_inversion_Wisting_cluster import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()
# mpirun -np 4 python main_cluster.py


def main():
    # rho_starting_model, rho_observed, recovered_model, recovered_data, data_observed = gravity_inversion_execution(
    # ) gravity_results_plotting_execution(rho_starting_model, rho_observed, recovered_model, recovered_data,
    # data_observed)
    dsw_starting_model, dsw_observed, data_observed, recovered_model, recovered_data, = \
        gravity_inversion_execution_Wisting()

    np.savetxt('dsw_starting_model_' + str(fig_num) + '.txt', dsw_starting_model, delimiter=' ', header='')
    np.savetxt('dsw_observed_' + str(fig_num) + '.txt', dsw_observed, delimiter=' ', header='')
    np.savetxt('data_observed_' + str(fig_num) + '.txt', data_observed, delimiter=' ', header='')
    np.savetxt('recovered_model_' + str(fig_num) + '.txt', recovered_model, delimiter=' ', header='')
    np.savetxt('recovered_data_' + str(fig_num) + '.txt', recovered_data, delimiter=' ', header='')


if __name__ == '__main__':
    main()
