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
        Vector<PatchInt> patches_int;
        Vector<PatchExt> patches_ext;

    public:
        // FV_Grid(Cells &&cells, Faces &&faces, Vector<PatchExt> &&PatchExtes, Vector<Patches> &&partition_PatchExtes);
        //  FV_Grid(Config &config, PrimalGrid &primal_grid);

        const Cells &get_cells() const { return cells; }
        const Faces &get_faces() const { return faces; }
        const Vector<PatchExt> &get_patches_ext() const { return patches_ext; }
        const Vector<PatchInt> &get_patches_int() const { return patches_int; }

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

        Index find_num_ghost_ext() const;
        Index find_num_ghost_int() const;

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