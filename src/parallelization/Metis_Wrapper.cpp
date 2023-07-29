#include "../../include/parallelization/Metis_Wrapper.hpp"

namespace NF_METIS
{
    Vector<Index> calc_element_partition(PrimalGrid &primal_grid,
                                         Index n_partitions)
    {
        static_assert(sizeof(idx_t) == sizeof(Index));

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
#ifndef NDEBUG
        options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; // Enable debug info
#endif

        Elements &elements = primal_grid.get_vol_elements();

        Vector<Index> &e_ptr = elements.get_e_ptr();
        Vector<Index> &e_ind = elements.get_e_ind();

        idx_t n_elements = static_cast<idx_t>(elements.size());
        idx_t n_nodes = primal_grid.get_nodes().size();

        idx_t n_common = 3; // Minimum 3 common node to get a face/edge for a 3D mesh (A triangle)
        idx_t objval = -1;
        Vector<Index> element_partition(n_elements);
        Vector<Index> node_partition(n_nodes); // Don't really need this one, but I get segfault if I input a nullptr instead of this
        idx_t n_part_idx_t = static_cast<idx_t>(n_partitions);

        int status = METIS_PartMeshDual(&n_elements,
                                        &n_nodes,
                                        (idx_t *)e_ptr.data(),
                                        (idx_t *)e_ind.data(),
                                        nullptr,
                                        nullptr,
                                        &n_common,
                                        &n_part_idx_t,
                                        nullptr,
                                        nullptr,
                                        &objval,
                                        (idx_t *)element_partition.data(),
                                        (idx_t *)node_partition.data());

        assert(metis_statuses.count((rstatus_et)status) == 1);
        if (status != 1)
            throw std::runtime_error("Error partitioning mesh using METIS. METIS_PartMeshDual return message: " +
                                     metis_statuses.at(static_cast<rstatus_et>(status)) + "\n");
        return element_partition;
    }
}