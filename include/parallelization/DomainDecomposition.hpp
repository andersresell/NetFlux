#pragma once
#include "../geometry/FV_Grid.hpp"
#include <metis.h>

namespace geometry
{
    class GridCreator
    {
        static constexpr int NOT_ASSIGNED{-1};

    public:
        static void create_partitioned_grids(Config &config,
                                             unique_ptr<PrimalGrid> &primal_grid,
                                             unique_ptr<FV_Grid> &FV_grid);

    private:
        static void construct_global_face_entities(const Elements &volume_elements_glob,
                                                   const Vector<ElementPatch> &element_patches,
                                                   map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                                   Elements &face_elements_glob,
                                                   map<FaceElement, InternalGhostData> &internal_boundary_faces,
                                                   const Vector<Index> &part);

        static void create_local_volume_entities(const Elements &volume_elements_glob,
                                                 const Vector<Vec3> &nodes_glob,
                                                 const Vector<ElementPatch> &element_patches_glob,
                                                 Index i_part,
                                                 const Vector<Index> &part,
                                                 Elements &volume_elements_loc,
                                                 Vector<Vec3> &nodes_loc,
                                                 Vector<ElementPatch> &element_patches_loc,
                                                 map<Index, Index> &nID_glob_to_loc,
                                                 map<Index, Index> &nID_loc_to_glob,
                                                 map<Index, Index> &eID_glob_to_loc);

        static void create_FV_grid_local(const map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                         const GridUtils &grid_utils,
                                         Index rank_loc,
                                         unique_ptr<FV_Grid> &FV_Grid_loc,
                                         PrimalGrid &primal_grid_loc,
                                         map<FaceElement, InternalGhostData> const &internal_boundary_faces_glob);
    };

    struct InternalGhostData
    {
        Index cID_a, cID_b;
        ShortIndex rank_a, rank_b;
    };

    class GridUtils
    {
        const Vector<Index> &part;
        const map<Index, Index> &eID_glob_to_loc;
        const map<Index, Index> &nID_glob_to_loc;
        const map<Index, Index> &nID_loc_to_glob;

    public:
        GridUtils(const Vector<Index> &part,
                  const map<Index, Index> &eID_glob_to_loc,
                  const map<Index, Index> &nID_glob_to_loc,
                  const map<Index, Index> &nID_loc_to_glob)
            : part{part}, eID_glob_to_loc{eID_glob_to_loc},
              nID_glob_to_loc{nID_glob_to_loc}, nID_loc_to_glob{nID_loc_to_glob}
        {
        }

        const map<Index, Index> get_map_nodeID_glob_to_loc() const
        {
            return nID_glob_to_loc;
        }
        const map<Index, Index> get_map_nodeID_loc_to_glob() const
        {
            return nID_loc_to_glob;
        }

        Index get_partID_from_global_eID(Index eID_glob) const
        {
            return part[eID_glob];
        }
        Index get_local_ghost_index(Index fID_glob, Index cID_i_glob) const
        {
            Index u = face_graph_glob.get_u(fID_glob);
            Index v = face_graph_glob.get_v(fID_glob);
            assert(u != v);

            Index ghostID_glob;
            if (cID_i_glob == u)
            {
                ghostID_glob = v;
            }
            else
            {
                assert(cID_i_glob == v);
                ghostID_glob = u;
            }
            assert(eID_glob_to_loc.count(ghostID_glob) == 1);
            return eID_glob_to_loc.at(ghostID_glob);
        }
        Index get_local_element_index(Index eID_glob) const
        {
            assert(eID_glob_to_loc.at(eID_glob) == 1);
            return eID_glob_to_loc.at(eID_glob);
        }
        Index get_local_node_index(Index nID_glob) const
        {
            assert(nID_glob_to_loc.at(nID_glob) == 1);
            return nID_glob_to_loc.at(nID_glob);
        }
    };

}

namespace NF_METIS
{
    using namespace geometry;

    const map<rstatus_et, string> metis_statuses = {{METIS_OK, "METIS_OK"},
                                                    {METIS_ERROR_INPUT, "METIS_ERROR_INPUT"},
                                                    {METIS_ERROR_MEMORY, "METIS_ERROR_MEMORY"},
                                                    {METIS_ERROR, "METIS_ERROR"}};

    Vector<Index> calc_element_partition(PrimalGrid &primal_grid,
                                         Index n_partitions);

}