#include "../../include/parallelization/Decomposition.hpp"

/*Assuming tetrahedron only for now*/
void METIS::create_mesh_partition(const Vector<TetConnect> &elements,
                                  Vector<idx_t> &element_partition,
                                  idx_t n_nodes,
                                  idx_t n_partitions)
{
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
#ifndef NDEBUG
    options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; // Enable debug info
#endif

    idx_t n_elements = elements.size();

    /*--------------------------------------------------------------------
    Represent the mesh connectivity on the format that METIS uses
    --------------------------------------------------------------------*/
    Vector<idx_t> eptr(n_elements + 1);
    Vector<idx_t> eind(n_elements * N_TET_NODES);

    eptr[0] = 0;
    for (Index i{0}; i < elements.size(); i++)
    {
        eptr[i + 1] = (i + 1) * N_TET_NODES;
        for (ShortIndex j{0}; j < N_TET_NODES; j++)
            eind[i * N_TET_NODES] = elements[i][j];
    }

    idx_t n_common = 3; // Minimum 3 common node to get a face/edge for a 3D mesh (A triangle)
    idx_t objval = -1;
    Vector<idx_t> node_partition(n_nodes); // Don't really need this one, but I get segfault if I input a nullptr instead of this
    int status = METIS_PartMeshDual(&n_elements,
                                    (idx_t *)&n_nodes,
                                    eptr.data(),
                                    eind.data(),
                                    nullptr,
                                    nullptr,
                                    &n_common,
                                    &n_partitions,
                                    nullptr,
                                    nullptr,
                                    &objval,
                                    element_partition.data(),
                                    node_partition.data());

    if (status != 1)
        throw std::runtime_error("Error partitioning mesh using METIS. METIS_PartMeshDual return message: " +
                                 metis_statuses.at(static_cast<rstatus_et>(status)) + "\n");
}