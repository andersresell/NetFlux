#pragma once
#include "../Grid.hpp"

namespace METIS
{
#include <metis.h>
    using namespace geom;

    const map<rstatus_et, string> metis_statuses = {{METIS_OK, "METIS_OK"},
                                                    {METIS_ERROR_INPUT, "METIS_ERROR_INPUT"},
                                                    {METIS_ERROR_MEMORY, "METIS_ERROR_MEMORY"},
                                                    {METIS_ERROR, "METIS_ERROR"}};

    void create_mesh_partition(const Vector<TetConnect> &elements,
                               Vector<idx_t> &element_partition,
                               idx_t n_nodes,
                               idx_t n_partitions);
}