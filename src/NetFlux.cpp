
#include "../include/parallelization/MPI.hpp"

#include "../include/Driver.hpp"

#include <mpi.h>

int main(int argc, char *argv[])
{

    int rank, size_world;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size_world);
    cout << "rank " << rank << ", size " << size_world << endl;

    try
    {
        if (argc != 2)
        {
            cout << "Specify config file to run NetFlux.\nUsage: $ ./NetFlux <path-to-sim-directory>\n";
            return 1;
        }
        string sim_dir_path = argv[1];
        cout << "\nProject folder: "
             << "'" << sim_dir_path << "'\n";
        Config config{sim_dir_path};
        Driver driver{config};
        driver.solve();
    }
    catch (std::exception &e)
    {
        std::cerr << "Exception caught:\n"
                  << e.what() << endl;
        return 1;
    }

    return 0;
}