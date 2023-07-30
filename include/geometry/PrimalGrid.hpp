#pragma once
#include "GeometryUtilities.hpp"
#include "../Config.hpp"

namespace geometry
{

    /*--------------------------------------------------------------------
    PrimalGrid stores reads the mesh file and stores nodes and
    connectivities of the FE type native/primal mesh, which is used to
    construct the face-based structure used by the FV solver.
    --------------------------------------------------------------------*/
    class PrimalGrid
    {
        friend class FV_Grid;
        friend class GridCreator;
        /*These are read from the input mesh*/
        Vector<Vec3> nodes;
        Elements vol_elements;
        Vector<ElementPatch> element_patches;
        Elements face_elements; /*Storing elements of all faces*/

        void read_mesh(const Config &config);
        void read_netflux_mesh(const Config &config);
        void read_su2_mesh(const Config &config);

    public:
        /*Reads mesh file*/
        PrimalGrid(const Config &config);
        /*Is constructed from existing mesh properties*/
        PrimalGrid(Vector<Vec3> &&nodes, Elements &&vol_elements, Vector<ElementPatch> &&element_patches);

        const Vector<Vec3> &get_nodes() const { return nodes; }
        Vector<Vec3> &get_nodes() { return nodes; }

        const Elements &get_vol_elements() const { return vol_elements; }
        Elements &get_vol_elements() { return vol_elements; }

        const Elements &get_face_elements() const { return face_elements; }
        Elements &get_face_elements() { return face_elements; }

        const Vector<ElementPatch> &get_element_patches() const { return element_patches; }
        Vector<ElementPatch> &get_element_patches() { return element_patches; }

        void print_grid() const;
        void partial_validity_check();

        Index find_num_ghost_external() const;

        template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar &vol_elements;
            ar &face_elements;
        }
    };

}