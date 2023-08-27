#pragma once
#include "PrimalGrid.hpp"
#include "../parallelization/Communicator.hpp"

namespace geometry
{

    /*--------------------------------------------------------------------
    FV_Grid contains mainly the face-based structure used by the FV solvers.
    --------------------------------------------------------------------*/
    class FV_Grid
    {
        friend class GridCreator;

        Cells cells;
        Faces faces;
        vector<PatchInterface> patches_interf;
        vector<PatchBoundary> patches_bound;

    public:
        FV_Grid() = default;
        // FV_Grid(Cells &&cells, Faces &&faces, vector<PatchBoundary> &&patches, vector<Patches> &&partition_patches);
        //  FV_Grid(Config &config, PrimalGrid &primal_grid);

        const Cells &get_cells() const { return cells; }
        const Faces &get_faces() const { return faces; }
        const vector<PatchBoundary> &get_patches_boundary() const { return patches_bound; }
        const vector<PatchInterface> &get_patches_interface() const { return patches_interf; }

        template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar &cells;
            ar &faces;
            ar &patches_interf;
            ar &patches_bound;
        }

    private:
        /*Creates computational grid from native mesh*/
        // void create_face_structure(Config &config, PrimalGrid &primal_grid);
        void resize_geometry_properties()
        {
            // cells.resize(cells.size());
            faces.resize_geometry_properties();
        }
        /*Calculates various geometrical properties to cells and faces based on the primal grid*/
        void calc_geometry_properties(const Config &config, const PrimalGrid &primal_grid, PartitionComm &part_comm);

        /*------Helper functions for creating face structure-------*/

        Index find_num_ghost_ext() const;
        Index find_num_ghost_part() const;

        void calc_face_properties(ElementType e_type,
                                  const Index *element,
                                  const vector<Vec3> &nodes,
                                  const Vec3 &cell_center_i,
                                  const Vec3 &cell_center_j,
                                  Vec3 &S_ij,
                                  Vec3 &centroid_to_face_i,
                                  Vec3 &centroid_to_face_j);

        void calc_ghost_centroid_bound(ElementType boundary_e_type,
                                       const Index *boundary_element,
                                       const vector<Vec3> &nodes,
                                       const Vec3 &centroid_i,
                                       Vec3 &centroid_ghost);
        void print_grid(const Config &config, const Elements &face_elements) const;
    };

};