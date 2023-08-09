#include "../../include/parallelization/MetisWrapper.hpp"

namespace NF_METIS
{
    vector<ShortIndex> calc_element_partition(PrimalGrid &primal_grid,
                                              Index n_partitions)
    {
        assert(NF_MPI::get_rank() == 0);
        static_assert(sizeof(idx_t) == sizeof(Index));

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
#ifndef NDEBUG
        options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; // Enable debug info
#endif

        Elements &elements = primal_grid.get_vol_elements();

        vector<Index> &e_ptr = elements.get_e_ptr();
        vector<Index> &e_ind = elements.get_e_ind();

        idx_t n_elements = static_cast<idx_t>(elements.size());
        idx_t n_nodes = primal_grid.get_nodes().size();

        idx_t n_common = 3; // Minimum 3 common node to get a face/edge for a 3D mesh (A triangle)
        idx_t objval;
        vector<idx_t> element_partition(n_elements, 0);
        vector<idx_t> node_partition(n_nodes); // Don't really need this one, but I get segfault if I input a nullptr instead of this
        idx_t n_part = static_cast<idx_t>(n_partitions);

        // Vector<idx_t> e_ptr_cp(e_ptr.size());
        // std::copy(e_ptr.begin(), e_ptr.end(), e_ptr_cp.begin());
        // Vector<idx_t> e_ind_cp(e_ind.size());
        // std::copy(e_ind.begin(), e_ind.end(), e_ind_cp.begin());

        // int status = METIS_PartMeshDual(&n_elements,
        //                                 &n_nodes,
        //                                 e_ptr_cp.data(),
        //                                 e_ind_cp.data(),
        //                                 NULL,
        //                                 NULL,
        //                                 &n_common,
        //                                 &n_part,
        //                                 NULL,
        //                                 options,
        //                                 &objval,
        //                                 element_partition.data(),
        //                                 node_partition.data());

        // idx_t eptr[] = {0, 3, 6, 9, 12};
        // idx_t eind[] = {0, 1, 2,
        //                 1, 3, 4,
        //                 1, 4, 2,
        //                 2, 4, 5};

        // idx_t ne = 4;
        // idx_t nn = 6;
        // idx_t ncommon = 2;
        // idx_t nparts = 4;
        // // idx_t objval;
        // vector<idx_t> epart(ne);
        // vector<idx_t> npart(nn);

        // idx_t options_[METIS_NOPTIONS];
        // METIS_SetDefaultOptions(options_);
        // options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; // Enable debug info

        // cerr << "starting metis\n";
        // int status = METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL,
        //                                 &ncommon, &nparts, NULL, options_, &objval, epart.data(), npart.data());

        if (n_part > 1)
        {
            int status = METIS_PartMeshDual(&n_elements,
                                            &n_nodes,
                                            (idx_t *)e_ptr.data(),
                                            (idx_t *)e_ind.data(),
                                            nullptr,
                                            nullptr,
                                            &n_common,
                                            &n_part,
                                            nullptr,
                                            nullptr,
                                            &objval,
                                            (idx_t *)element_partition.data(),
                                            (idx_t *)node_partition.data());

            assert(metis_statuses.count((rstatus_et)status) == 1);
            if (status != 1)
                throw std::runtime_error("Error partitioning mesh using METIS. METIS_PartMeshDual return message: " +
                                         metis_statuses.at(static_cast<rstatus_et>(status)) + "\n");
            vector<ShortIndex> e_part_short;
            e_part_short.reserve(element_partition.size());
            for (Index rank : element_partition)
                e_part_short.emplace_back(rank);

            return e_part_short;
        }
        else
        {
            assert(n_part == 1);
            return vector<ShortIndex>(element_partition.size(), 0);
        }
    }
}