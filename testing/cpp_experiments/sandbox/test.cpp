#include <iostream>
#include <cassert>
#include <limits>
using namespace std;

#include "../../../include/SolverData.hpp"
#include "../../../include/containers/DynamicContainer.hpp"

#include "test.hpp"

int main(int argc, char *argv[])
{

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the rank and size of the MPI communicator
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Print "Hello, World!" along with the rank of the process
    std::cout << "Hello, World! I am process " << rank << " out of " << size << " processes." << std::endl;

    // Finalize MPI
    MPI_Finalize();

    return 0;
}