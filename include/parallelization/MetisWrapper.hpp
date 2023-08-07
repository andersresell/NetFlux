#pragma once
#include "../geometry/FV_Grid.hpp"
#include <metis.h>
namespace NF_METIS
{
    using namespace geometry;

    const map<rstatus_et, string> metis_statuses = {{METIS_OK, "METIS_OK"},
                                                    {METIS_ERROR_INPUT, "METIS_ERROR_INPUT"},
                                                    {METIS_ERROR_MEMORY, "METIS_ERROR_MEMORY"},
                                                    {METIS_ERROR, "METIS_ERROR"}};

    vector<ShortIndex> calc_element_partition(PrimalGrid &primal_grid,
                                              Index n_partitions);

}