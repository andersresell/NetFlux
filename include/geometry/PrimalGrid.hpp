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
        vector<Vec3> nodes;
        Elements vol_elements;
        vector<ElementPatch> element_patches;
        Elements face_elements; /*Storing elements of all faces*/
        // const Index eID_glob_first{0}; /*global element index of first element*/

        void read_mesh(const Config &config);
        void read_netflux_mesh(const Config &config);
        void read_su2_mesh(const Config &config);

    public:
        /*Reads mesh file*/
        PrimalGrid() = default;

        PrimalGrid(const Config &config);
        /*Is constructed from existing mesh properties*/
        PrimalGrid(vector<Vec3> &&nodes, Elements &&vol_elements, vector<ElementPatch> &&element_patches);

        const vector<Vec3> &get_nodes() const { return nodes; }
        vector<Vec3> &get_nodes() { return nodes; }

        const Elements &get_vol_elements() const { return vol_elements; }
        Elements &get_vol_elements() { return vol_elements; }

        const Elements &get_face_elements() const { return face_elements; }
        Elements &get_face_elements() { return face_elements; }

        const vector<ElementPatch> &get_element_patches() const { return element_patches; }
        vector<ElementPatch> &get_element_patches() { return element_patches; }

        // Index get_eID_global_first() const { return eID_glob_first; }

        void print_grid() const;
        void partial_validity_check();

        Index find_num_ghost_external() const;

        template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar &nodes;
            ar &vol_elements;
            ar &face_elements;
            ar &element_patches;
            ar &face_elements;
        }
    };

}