#pragma once
#include "MPI_Wrapper.hpp"
#include "../geometry/GeometryUtilities.hpp"
#include "../containers/DynamicContainer.hpp"

/*Used to send or receive data across a patch. Contains enough storage
to send all the fields store all the */

template <ShortIndex N_COLS>
using Field = DynamicContainer3D<Scalar, N_COLS>;

class InterfaceComm
{
    // const map<FieldType, ShortIndex> field_type_dimension = {{FieldType::Primvars, 1},
    //                                                          {FieldType::PrimvarsGrad, N_DIM},
    //                                                          {FieldType::Consvars, 1}};
    // const ShortIndex N_EQS;
    vector<Scalar> sendbuf;
    vector<Scalar> recvbuf;
    const geometry::PatchInterface &patch_part;
    const geometry::Faces &faces;
    Index sendptr, recvptr;
    // //const Index size_vecfield, size_gradfield;
    // const ShortIndex n_vecfields_max, n_gradfields_max;

public:
    InterfaceComm(const geometry::PatchInterface &patch_part,
                  const geometry::Faces &faces);

    void send_receive_fields(MPI_Request &send_req, MPI_Request &recv_req);

    Index n_ghost() const { return patch_part.N_FACES; }

    void set_max_size_cell(ShortIndex max_scalars_per_cell)
    {
        sendbuf.resize(n_ghost() * max_scalars_per_cell);
        recvbuf.resize(n_ghost() * max_scalars_per_cell);
    }
    ShortIndex rank_neigbour() const { return patch_part.rank_neighbour; }

    void clear()
    {
        assert(sendptr == recvptr); /*These should be equal with correct use*/
        sendptr = 0;
        recvptr = 0;
    }

    bool is_clear() { return sendptr == 0 && recvptr == 0; }

    template <ShortIndex N_COLS>
    void pack_field(const Field<N_COLS> &sendfield)
    {
        static_assert(N_COLS == 1 || N_COLS == N_DIM);
        for (Index ij{patch_part.FIRST_FACE}; ij < patch_part.FIRST_FACE + patch_part.N_FACES; ij++)
        {
            Index i_domain = faces.get_cell_i(ij);
            for (ShortIndex neq{0}; neq < sendfield.rows(); neq++)
                for (ShortIndex idim{0}; idim < N_COLS; idim++)
                    sendbuf[sendptr++] = sendfield(i_domain, neq, idim);
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
            for (ShortIndex neq{0}; neq < recvfield.rows(); neq++)
                for (ShortIndex idim{0}; idim < N_COLS; idim++)
                    recvfield(j_ghost, neq, idim) = recvbuf[recvptr++];
        }
        assert(recvptr < recvbuf.size());
    }

    void pack_Vec3_field(const vector<Vec3> &sendfield)
    {
        for (Index ij{patch_part.FIRST_FACE}; ij < patch_part.FIRST_FACE + patch_part.N_FACES; ij++)
        {
            Index i_domain = faces.get_cell_i(ij);
            for (ShortIndex idim{0}; idim < N_DIM; idim++)
                sendbuf[sendptr++] = sendfield[i_domain][idim];
        }
        assert(sendptr <= sendbuf.size());
    };

    void unpack_Vec3_field(vector<Vec3> &recvfield)
    {
        for (Index ij{patch_part.FIRST_FACE}; ij < patch_part.FIRST_FACE + patch_part.N_FACES; ij++)
        {
            Index j_ghost = faces.get_cell_j(ij);
            for (ShortIndex idim{0}; idim < N_DIM; idim++)
                recvfield[j_ghost][idim] = recvbuf[recvptr++];
        }
        assert(recvptr <= recvbuf.size());
    };
};

/*Handles sending and receiving of fields to ghost cells at the interfaces
between local domains/partitions.*/
class PartitionComm
{
    // const Index N_CELLS_TOT;

    vector<unique_ptr<InterfaceComm>> interf_comms;
    // const ShortIndex n_vecfields_max, n_gradfields_max;
    ShortIndex max_scalars_per_cell_;

public:
    PartitionComm(const Config &config,
                  const geometry::Faces &faces,
                  const vector<geometry::PatchInterface> &patches_interf);

    void set_max_size_cell(ShortIndex max_scalars_per_cell)
    {
        max_scalars_per_cell_ = max_scalars_per_cell;
        for (const auto &ic : interf_comms)
            ic->set_max_size_cell(max_scalars_per_cell);
    }

    ShortIndex get_max_size_cell() const
    {
        return max_scalars_per_cell_;
    }

    void communicate_ghost_fields();
    void clear();
    bool is_clear();

    // Index get_n_vecfields_max() const { return n_vecfields_max; }
    // Index get_n_gradfields_max() const { return n_gradfields_max; }

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

    void pack_Vec3_field(const vector<Vec3> &sendfield)
    {
        for (auto &interf_comm : interf_comms)
            interf_comm->pack_Vec3_field(sendfield);
    };

    void unpack_Vec3_field(vector<Vec3> &recvfield)
    {
        for (auto &interf_comm : interf_comms)
            interf_comm->unpack_Vec3_field(recvfield);
    };

    void communicate_interface_ghost_centroids(vector<Vec3> &centroids);

private:
    ShortIndex num_patches() const { return interf_comms.size(); }
};
