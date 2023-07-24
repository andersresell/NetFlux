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
        /*These are read from the input mesh*/
        Vector<Vec3> nodes;
        Elements volume_elements;
        Vector<ElementPatch> element_patches;

        Elements face_elements; /*Storing elements of all faces*/

        void read_mesh(const Config &config);
        void read_netflux_mesh(const Config &config);
        void read_su2_mesh(const Config &config);

    public:
        PrimalGrid(const Config &config);
        const Vector<Vec3> &get_nodes() const { return nodes; }
        const Elements &get_volume_elements() const { return volume_elements; }
        const Elements &get_face_elements() const { return face_elements; }
        const Vector<ElementPatch> &get_element_patches() const { return element_patches; }
        Elements &get_face_elements() { return face_elements; }
        void print_grid() const;
    };

}