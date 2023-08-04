#pragma once
#include "MPI_Wrapper.hpp"
#include "../geometry/FV_Grid.hpp"

/*Used to send or receive data across a patch. Contains enough storage
to send all the fields store all the */

template <ShortIndex N_COLS>
using Field = DynamicContainer3D<Scalar, N_COLS>;

class PartComm
{
    // const map<FieldType, ShortIndex> field_type_dimension = {{FieldType::Primvars, 1},
    //                                                          {FieldType::PrimvarsGrad, N_DIM},
    //                                                          {FieldType::Consvars, 1}};
    const ShortIndex N_EQS;
    Vector<Scalar> sendbuf;
    Vector<Scalar> recvbuf;
    const geometry::PatchPart &patch_part;
    const geometry::Faces &faces;
    Index sendptr, recvptr;
    const Index size_vecfield, size_gradfield;
    const ShortIndex n_vecfields_max, n_gradfields_max;
    bool comm_complete_{false};

public:
    PartComm(const ShortIndex n_vecfields_max,
             const ShortIndex n_gradfields_max,
             const ShortIndex N_EQS,
             const geometry::PatchPart &patch_part,
             const geometry::Faces &faces);

    void send_and_receive_fields();
    bool comm_complete() const { return comm_complete_; }

    Index n_ghost() const { return patch_part.N_FACES; }
    ShortIndex rank_neigbour() const { return patch_part.rank_neighbour; }

    void clear()
    {
        sendptr = 0;
        recvptr = 0;
        comm_complete_ = false;
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
