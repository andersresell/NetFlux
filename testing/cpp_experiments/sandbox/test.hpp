#pragma once
#include "../../include/Includes.hpp"

struct Faces
{
    friend class Grid;

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
    friend class Grid;
    Vector<Scalar> cell_volumes;
    Vector<Vec3> centroids;

    void reserve(Index size);
    void add_empty();

public:
    Index size() const { return cell_volumes.size(); }
    Scalar get_cell_volume(Index cell_index) const { return cell_volumes[cell_index]; }
    const Vec3 &get_centroid(Index cell_index) const { return centroids[cell_index]; }
};
