#pragma once
#include "MPI_Wrapper.hpp"
#include "../geometry/FV_Grid.hpp"

/*Used to send or receive data across a patch. Contains enough storage
to send all the fields store all the */

template <ShortIndex N_COLS>
using Field = DynamicContainer3D<Scalar, N_COLS>;

class InterfaceComm
{
    // const map<FieldType, ShortIndex> field_type_dimension = {{FieldType::Primvars, 1},
    //                                                          {FieldType::PrimvarsGrad, N_DIM},
    //                                                          {FieldType::Consvars, 1}};
    const ShortIndex N_EQS;
    Vector<Scalar> sendbuf;
    Vector<Scalar> recvbuf;
    const geometry::PatchInterface &patch_part;
    const geometry::Faces &faces;
    Index sendptr, recvptr;
    // //const Index size_vecfield, size_gradfield;
    // const ShortIndex n_vecfields_max, n_gradfields_max;

public:
    InterfaceComm(ShortIndex n_vecfields_max,
                  ShortIndex n_gradfields_max,
                  const ShortIndex N_EQS,
                  const geometry::PatchInterface &patch_part,
                  const geometry::Faces &faces);

    void send_receive_fields(MPI_Request &send_req, MPI_Request &recv_req);

    Index n_ghost() const { return patch_part.N_FACES; }
    ShortIndex rank_neigbour() const { return patch_part.rank_neighbour; }

    void clear()
    {
        assert(sendptr == recvptr); /*These should be equal with correct use*/
        sendptr = 0;
        recvptr = 0;
    }

    template <ShortIndex N_COLS>
    void pack_field(const Field<N_COLS> &sendfield)
    {
        static_assert(N_COLS == 1 || N_COLS == N_DIM);
        for (Index ij{patch_part.FIRST_FACE}; ij < patch_part.FIRST_FACE + patch_part.N_FACES; ij++)
        {
            Index i_domain = faces.get_cell_i(ij);
            for (ShortIndex neq{0}; neq < N_EQS; neq++)
                for (ShortIndex ndim{0}; ndim < N_COLS; ndim++)
                    sendbuf[sendptr++] = sendfield(i_domain, neq, ndim);
        }
        assert(sendptr < sendbuf.size());
    }

    template <ShortIndex N_COLS>
    void unpack_field(Field<N_COLS> &recvfield)
    {
        static_assert(N_COLS == 1 || N_COLS == N_DIM);
        for (Index ij{patch_part.FIRST_FACE}; ij < patch_part.FIRST_FACE + patch_part.N_FACES; ij++)
        {
            Index j_ghost = faces.get_cell_j(ij);
            for (ShortIndex neq{0}; neq < N_EQS; neq++)
                for (ShortIndex ndim{0}; ndim < N_COLS; ndim++)
                    recvfield(j_ghost, neq, n_dim) = recvbuf[recvptr++];
        }
        assert(recvptr < recvbuf.size());
    }
};

/*Handles sending and receiving of fields to ghost cells at the interfaces
between local domains/partitions.*/
class PartitionComm
{
    Vector<unique_ptr<InterfaceComm>> interf_comms;

    const ShortIndex n_vecfields_max, n_gradfields_max;

public:
    PartitionComm(const ShortIndex n_vecfields_max,
                  const ShortIndex n_gradfields_max,
                  const ShortIndex N_EQS,
                  const geometry::FV_Grid &FV_grid);

    void communicate_ghost_fields();
    void clear();

    Index get_n_vecfields_max() const { return n_vecfields_max; }
    Index get_n_gradfields_max() const { return n_gradfields_max; }

    template <ShortIndex N_COLS>
    void pack_field(const Field<N_COLS> &sendfield)
    {
        for (auto &interf_comm : interf_comms)
            interf_comm->pack_field(sendfield);
    }

    template <ShortIndex N_COLS>
    void unpack_field(Field<N_COLS> &recvfield)
    {
        for (auto &interf_comm : interf_comms)
            interf_comm->unpack_field(recvfield);
    }

private:
    ShortIndex num_patches() const { return interf_comms.size(); }
};