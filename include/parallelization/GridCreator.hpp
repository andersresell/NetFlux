#pragma once
#include "../geometry/FV_Grid.hpp"
#include "Metis_Wrapper.hpp"
#include "Serialization.hpp"

namespace geometry
{
    class GridCreator
    {
        static constexpr int CELL_NOT_ASSIGNED{-1};

    public:
        static void create_partitioned_grids(Config &config,
                                             unique_ptr<PrimalGrid> &primal_grid,
                                             unique_ptr<FV_Grid> &FV_grid);

    private:
        static void reorder_global_grid(const Vector<Index> &part,
                                        PrimalGrid &primal_grid,
                                        Vector<pair<Index, Index>> &part_to_element_range,
                                        Vector<Index> &eID_glob_to_loc);

        static void create_global_face_entities(const Elements &vol_elements_glob,
                                                const Vector<ElementPatch> &element_patches,
                                                map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                                Elements &face_elements_glob,
                                                map<FaceElement, GhostDataPartition> &internal_boundary_faces,
                                                const Vector<Index> &part);

        static void create_primal_grid_local(const Elements &vol_elements_glob,
                                             const Vector<Vec3> &nodes_glob,
                                             const Vector<ElementPatch> &element_patches_glob,
                                             Index rank_loc,
                                             const Vector<pair<Index, Index>> &part_to_element_range,
                                             const Vector<Index> &part,
                                             unique_ptr<PrimalGrid> &primal_grid_loc,
                                             Vector<Index> &eID_glob_to_loc,
                                             map<Index, Index> &nID_glob_to_loc,
                                             map<Index, Index> &nID_loc_to_glob);

        static void create_FV_grid_local(const Config &config,
                                         const map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                         const PartitionUtils &utils,
                                         Index rank_loc,
                                         unique_ptr<FV_Grid> &FV_grid_loc,
                                         PrimalGrid &primal_grid_loc,
                                         map<FaceElement, GhostDataPartition> const &internal_boundary_faces_glob);

        /*Sorts faces and face_elements patch-wise after the cell indices of the faces.*/
        static void reorder_faces_enitities(Index num_interior_faces,
                                            const Vector<Patch> &patches,
                                            Faces &faces,
                                            Elements &face_elements);

        static void send_recv_grids(Config &config,
                                    Vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                    Vector<unique_ptr<FV_Grid>> &FV_grids_loc,
                                    unique_ptr<PrimalGrid> &primal_grid,
                                    unique_ptr<FV_Grid> &FV_grid);
        static void set_config_grid_data_local(Config &config,
                                               unique_ptr<PrimalGrid> &primal_grid,
                                               unique_ptr<FV_Grid> &FV_grid);
    };

    struct GhostDataPartition
    {
        Index cID_a, cID_b;
        ShortIndex rank_a, rank_b;
    };

    class PartitionUtils
    {
        const Vector<Index> &part;
        const Vector<pair<Index, Index>> &part_to_e_range;
        const Vector<Index> &eID_glob_to_loc;
        const map<Index, Index> &nID_glob_to_loc;
        const map<Index, Index> &nID_loc_to_glob;

    public:
        PartitionUtils(const Vector<Index> &part,
                       const Vector<pair<Index, Index>> &part_to_e_range,
                       const Vector<Index> &eID_glob_to_loc,
                       const map<Index, Index> &nID_glob_to_loc,
                       const map<Index, Index> &nID_loc_to_glob)
            : part{part}, part_to_e_range{part_to_e_range}, eID_glob_to_loc{eID_glob_to_loc},
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

        Index get_local_element_index(Index eID_glob) const { return eID_glob_to_loc[eID_glob]; }

        Index get_local_node_index(Index nID_glob) const
        {
            {
                assert(nID_glob_to_loc.at(nID_glob) == 1);
                return nID_glob_to_loc.at(nID_glob);
            }
        }
    };

}
