#pragma once
#include "../Includes.hpp"
#include <mpi.h>

class NF_MPI
{
    inline static int rank, size;

public:
    static void Init(int *argc, char ***argv)
    {
        MPI_Init(argc, argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }
    static ShortIndex get_rank() { return rank; }

    static ShortIndex get_size() { return size; }
};
