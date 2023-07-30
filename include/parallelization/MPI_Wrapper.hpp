#pragma once
#include "../Includes.hpp"
#include <mpi.h>

/*--------------------------------------------------------------------
Wrapping the necessary functionality of the C MPI API, removing
unneeded options
--------------------------------------------------------------------*/
class NF_MPI
{
    inline static int rank, size;

    template <typename T>
    static MPI_Datatype get_MPI_Datatype()
    {
        using std::is_same;
        if (is_same<T, double>::value)
            return MPI_DOUBLE;
        else if (is_same<T, uint32_t>::value)
            return MPI_UINT32_T;
        else if (is_same<T, char>::value)
            return MPI_BYTE;
        else
            static_assert(false);
    }

public:
    static void Init(int *argc, char ***argv)
    {
        MPI_Init(argc, argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }
    static ShortIndex get_rank() { return rank; }

    static ShortIndex get_size() { return size; }

    static void Finalize() { MPI_Finalize(); }

    static void Barrier() { MPI_Barrier(MPI_COMM_WORLD); }

    template <typename T>
    static void Send(const T *data, Index count, ShortIndex dest_rank, ShortIndex tag = 0)
    {
        MPI_Send(data, sizeof(T) * count, get_MPI_Datatype<T>(), dest_rank, tag, MPI_COMM_WORLD);
    }
    template <typename T>
    static void Recv(const T *data, Index count, ShortIndex source_rank, ShortIndex tag = 0)
    {
        MPI_Recv(data, count, get_MPI_Datatype<T>(), source_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
    }
};
