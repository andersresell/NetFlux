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
    static void Send(const T *sendbuf, Index count, ShortIndex dest_rank, ShortIndex tag = 0)
    {
        MPI_Send(sendbuf, sizeof(T) * count, get_MPI_Datatype<T>(), dest_rank, tag, MPI_COMM_WORLD);
    }
    template <typename T>
    static void Recv(const T *recvbuf, Index count, ShortIndex source_rank, ShortIndex tag = 0)
    {
        MPI_Recv(recvbuf, count, get_MPI_Datatype<T>(), source_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
    }

    template <typename T>
    static void Bcast(const T *sendbuf, Index count, ShortIndex source_rank)
    {
        MPI_Bcast(sendbuf, count, get_MPI_Datatype<T>(), source_rank, MPI_COMM_WORLD);
    }

    template <typename T>
    static void Gather(const T *sendbuf, Index sendcount, T *recvbuf, Index recvcount, ShortIndex dest_rank)
    {
        MPI_Gather(sendbuf, sendcount, get_MPI_Datatype<T>(), recvbuf, recvcount, get_MPI_Datatype<T>(), dest_rank, MPI_COMM_WORLD);
    }
};
