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

#include <iostream>
#include <array>

    // Raw array
    int rawArray[] = {1, 2, 3, 4, 5};

    // Length of the raw array
    size_t length = sizeof(rawArray) / sizeof(rawArray[0]);

    // Initialize the std::array from the raw pointer
    std::array<int, 5> myArray{rawArray, rawArray + length};

    // Now, myArray contains the elements from the raw array
    // You can use myArray just like any other std::array

    // Print the elements of the std::array
    for (int num : myArray)
    {
        std::cout << num << " ";
    }
    std::cout << std::endl;

    return 0;
}
