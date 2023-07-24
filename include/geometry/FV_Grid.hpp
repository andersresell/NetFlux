#pragma once
#include "PrimalGrid.hpp"

namespace geometry
{
    class Faces;
    class Cells;
    class Patches;

    /*--------------------------------------------------------------------
    FV_Grid contains mainly the face-based structure used by the FV solvers.
    --------------------------------------------------------------------*/
    class FV_Grid
    {
        Cells cells;
        Faces faces;
        Vector<Patch> patches;

    public:
        FV_Grid(Config &config, PrimalGrid &primal_grid);

        // void print_grid(const Config &config) const;
        void print_native_mesh() const;

        const Cells &get_cells() const { return cells; }
        const Faces &get_faces() const { return faces; }
        const Vector<Patch> &get_patches() const { return patches; }

    private:
        /*Creates computational grid from native mesh*/
        void create_face_structure(Config &config, PrimalGrid &primal_grid);

        /*Calculates various geometrical properties to cells and faces based on the primal grid*/
        void calc_geometry_properties(const Config &config, const PrimalGrid &primal_grid);

        /*------Helper functions for creating face structure-------*/

        // Reorders the (for now interior) faces in an optimal fashion based on the face indices
        void reorder_faces(const Config &config);

        /*Used to find the number of ghost cells before this value is set in Config object*/
        Index find_N_GHOST_cells(const Vector<ElementPatch> &element_patches);

        void calc_face_properties(ElementType e_type,
                                  const Index *element,
                                  const Vector<Vec3> &nodes,
                                  const Vec3 &cell_center_i,
                                  const Vec3 &cell_center_j,
                                  Vec3 &S_ij,
                                  Vec3 &centroid_to_face_i,
                                  Vec3 &centroid_to_face_j);

        void calc_ghost_centroid(ElementType face_e_type,
                                 const Index *surface_element,
                                 const Vector<Vec3> &nodes,
                                 const Vec3 &centroid_i,
                                 Vec3 &centroid_ghost);
    };

};