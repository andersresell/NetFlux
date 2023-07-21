#pragma once

#include "Includes.hpp"
#include "Utilities.hpp"
#include "Config.hpp"

namespace geom
{

    class Grid
    {

        // Native mesh

        // Computational grid

    public:
        Grid(Config &config);

        // void print_grid(const Config &config) const;
        void print_native_mesh() const;

        const Vector<Vec3> &get_nodes() const { return nodes; }
        const Vector<TetConnect> &get_tet_connect() const { return tet_connect; }
        const Cells &get_cells() const { return cells; }
        const Faces &get_faces() const { return faces; }
        const Vector<Patch> &get_patches() const { return patches; }

    private:
        /*Read mesh file. This populates the:
        - nodes -> nodes of native mesh
        - tet_connect -> tetrahedral elements connectivity
        - tri_patch_connect -> patches of boundary triangles and respective BC types
        */
        void read_mesh(const Config &config);
        void read_netflux_mesh(const Config &config);
        void read_su2_mesh(const Config &config);

        /*Creates computational grid from native mesh*/
        void create_grid(Config &config);

        /*------Helper functions for creating grid-------*/

        // Reorders the (for now interior) faces in an optimal fashion based on the face indices
        void reorder_faces(const Config &config);

        /*Assigns cell centers, boundary normals, etc*/
        void assign_geometry_properties(const Config &config, const Vector<Triangle> &face_triangles);

        /*Finds the neighbouring cell j of face ij of cell i. If no neigbour exist (boundary), it returns false*/
        std::pair<Index, bool> find_neigbouring_cell(Index i,
                                                     TriConnect face_ij,
                                                     const Vector<TetConnect> &tet_connect) const;

        Tetrahedron tet_from_connect(const TetConnect &tc) const;
        Triangle tri_from_connect(const TriConnect &tc) const;

        /*Used to find the number of ghost cells before this value is set in Config object*/
        Index find_N_GHOST_cells();

        /*Calling shrink_to_fit on the different members after grid construction to reduce allocated memory.
        (This function is probably redundant, since most vector sizes are set and not grown dynamically)*/
        void shrink_vectors();
    };

}
