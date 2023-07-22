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
        FV_Grid(Config &config, const PrimalGrid &primal_grid);

        // void print_grid(const Config &config) const;
        void print_native_mesh() const;

        const Cells &get_cells() const { return cells; }
        const Faces &get_faces() const { return faces; }
        const Vector<Patch> &get_patches() const { return patches; }

    private:
        /*Creates computational grid from native mesh*/
        void create_face_structure(Config &config, PrimalGrid &primal_grid);

        /*------Helper functions for creating grid-------*/

        // Reorders the (for now interior) faces in an optimal fashion based on the face indices
        void reorder_faces(const Config &config);

        /*Assigns cell centers, boundary normals, etc*/
        void assign_geometry_properties(const Config &config,
                                        const Elements &elements,
                                        const Vector<Vec3> &nodes,
                                        const Vector<Triangle> &face_triangles);

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

    /*A structure of arrays (SoA) containing the faces and their required properties*/
    class Faces
    {
        friend class FV_Grid;

        struct CellPair
        {
            CellPair(Index i, Index j) : i{i}, j{j} {} // For some reason I had to define the constructor to make emplace_back work

            Index i, j;
            bool operator<(CellPair rhs) const
            {
                if (i != rhs.i)
                    return i < rhs.i;

                assert(j != rhs.j); // This would imply that cell indices are identical

                return j < rhs.j;
            }
        };

        Vector<CellPair> cell_indices;
        Vector<Vec3> normal_areas;
        Vector<Vec3> centroid_to_face_i;
        Vector<Vec3> centroid_to_face_j;

        void reserve(Index size);
        void resize_geometry_properties();
        void sort(Index first, Index last);

    public:
        Index size() const { return cell_indices.size(); }

        Index get_cell_i(Index face_index) const { return cell_indices[face_index].i; }
        Index get_cell_j(Index face_index) const { return cell_indices[face_index].j; }
        const Vec3 &get_normal_area(Index face_index) const { return normal_areas[face_index]; }
        const Vec3 &get_centroid_to_face_i(Index face_index) const { return centroid_to_face_i[face_index]; }
        const Vec3 &get_centroid_to_face_j(Index face_index) const { return centroid_to_face_j[face_index]; }
    };

    /*A structure of arrays (SoA) containing the cells and their required properties*/
    class Cells
    {
        friend class FV_Grid;
        Vector<Scalar> cell_volumes;
        Vector<Vec3> centroids;

        void resize(Index size);
        void add_empty();

    public:
        Index size() const { return cell_volumes.size(); }
        Scalar get_cell_volume(Index cell_index) const { return cell_volumes[cell_index]; }
        const Vec3 &get_centroid(Index cell_index) const { return centroids[cell_index]; }
    };

    struct Patch
    {
        BoundaryType boundary_type;
        // Vector<Index> boundary_face_indices;
        Index N_FACES;
        Index FIRST_FACE;
    };

    void set_cell_properties(Index cell_index,
                             ElementType type,
                             const Index *element,
                             const Vector<Vec3> &nodes,
                             Vector<Scalar> &cell_volumes,
                             Vector<Vec3> &centroids);
};