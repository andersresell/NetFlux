#pragma once
#include "PrimalGrid.hpp"

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
        Vector<Patch> patches;
        Vector<PartitionPatch> partition_patches;

    public:
        FV_Grid(Cells &&cells, Faces &&faces, Vector<Patch> &&patches, Vector<PartitionPatch> &&partition_patches);
        // FV_Grid(Config &config, PrimalGrid &primal_grid);

        const Cells &get_cells() const { return cells; }
        const Faces &get_faces() const { return faces; }
        const Vector<Patch> &get_patches() const { return patches; }
        const Vector<PartitionPatch> &get_partition_patches() const { return partition_patches; }

        template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar &cells;
            ar &faces;
        }

    private:
        /*Creates computational grid from native mesh*/
        // void create_face_structure(Config &config, PrimalGrid &primal_grid);

        /*Calculates various geometrical properties to cells and faces based on the primal grid*/
        void calc_geometry_properties(const Config &config, const PrimalGrid &primal_grid);

        /*------Helper functions for creating face structure-------*/

        Index find_num_ghost_external() const;
        Index find_num_ghost_tot() const;

        void calc_face_properties(ElementType e_type,
                                  const Index *element,
                                  const Vector<Vec3> &nodes,
                                  const Vec3 &cell_center_i,
                                  const Vec3 &cell_center_j,
                                  Vec3 &S_ij,
                                  Vec3 &centroid_to_face_i,
                                  Vec3 &centroid_to_face_j);

        void calc_ghost_centroid(ElementType boundary_e_type,
                                 const Index *boundary_element,
                                 const Vector<Vec3> &nodes,
                                 const Vec3 &centroid_i,
                                 Vec3 &centroid_ghost);
        void print_grid(const Config &config, const Elements &face_elements) const;
    };

};