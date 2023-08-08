#include <iostream>
#include <vector>
#include <mpi.h>
#include <unistd.h>

using namespace std;

int main(int argc, char **argv)
{
    {
        int i = 0;
        while (0 == i)
            sleep(5);
    }
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "hey rank " << rank << endl;

    MPI_Finalize();
}
