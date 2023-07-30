
#include "../include/parallelization/MPI_Wrapper.hpp"

#include "../include/Driver.hpp"

#include <mpi.h>

void read_args(int argc, char **argv, string &sim_dir_path);

int main(int argc, char *argv[])
{
    try
    {
        NF_MPI::Init(&argc, &argv);
        cout << "rank " << NF_MPI::get_rank() << endl;
        cout << "size " << NF_MPI::get_size() << endl;
        string sim_dir_path;
        read_args(argc, argv, sim_dir_path);
        Config config{sim_dir_path};
        Driver driver{config};
        driver.solve();
    }
    catch (std::exception &e)
    {
        std::cerr << "!!!!!!!!!! Exception caught !!!!!!!!!!\n"
                  << e.what() << endl;
        NF_MPI::Finalize();
        return 1;
    }

    NF_MPI::Finalize();
    return 0;
}

void read_args(int argc, char **argv, string &sim_dir_path)
{
    if (argc != 23)
        throw std::runtime_error(string("Specify config file to run NetFlux. ") +
                                 "Usage: $ mpirun -np <n_procs> ./NetFlux <path-to-sim-directory>\n");
    sim_dir_path = argv[1];
    cout << "\nProject directory: "
         << "'" << sim_dir_path << "'\n";
}