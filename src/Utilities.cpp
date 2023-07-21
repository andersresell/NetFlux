#include "../include/Utilities.hpp"

namespace geom
{

    void Faces::reserve(Index size)
    {
        cell_indices.reserve(size);
        normal_areas.reserve(size);
        centroid_to_face_i.reserve(size);
        centroid_to_face_j.reserve(size);
    }

    void Faces::resize_geometry_properties()
    {
        normal_areas.resize(cell_indices.size());
        centroid_to_face_i.resize(cell_indices.size());
        centroid_to_face_j.resize(cell_indices.size());
    }

    /*Sorting all the Vectors from indices begin to end based on the cell_indices*/
    void Faces::sort(Index begin, Index end)
    {
        assert(begin < end && end <= cell_indices.size());

        Vector<Index> indices(end - begin);
        for (Index i{begin}; i < end; i++)
            indices[i - begin] = i;

        std::sort(indices.begin(), indices.end(), [this](Index a, Index b)
                  { return cell_indices[a] < cell_indices[b]; });
        for (Index i{begin}; i < end; i++)
        {
            std::swap(cell_indices[i], cell_indices[indices[i - begin]]);
            std::swap(normal_areas[i], normal_areas[indices[i - begin]]);
            std::swap(centroid_to_face_i[i], centroid_to_face_i[indices[i - begin]]);
            std::swap(centroid_to_face_j[i], centroid_to_face_j[indices[i - begin]]);
        }
    }

    void Cells::reserve(Index size)
    {
        cell_volumes.reserve(size);
        centroids.reserve(size);
    }

    void Cells::add_empty()
    {
        cell_volumes.emplace_back();
        centroids.emplace_back();
    }

    TriConnect tet_face_connectivity(TetConnect tc, ShortIndex face_k)
    {
        assert(face_k >= 0 && face_k < N_TET_FACES);
        if (face_k == 0)
            return {tc.a(), tc.b(), tc.c()};
        else if (face_k == 1)
            return {tc.a(), tc.b(), tc.d()};
        else if (face_k == 2)
            return {tc.a(), tc.c(), tc.d()};
        else
            return {tc.b(), tc.c(), tc.d()};
    }

    Vec3 Triangle::calc_area_normal() const
    {
        assert(nodes.size() == 3);
        // Using that the cross product of two vectors is the area of a parallelogram (with correct normal), so the half is the area of a triangle
        Vec3 ab = nodes.at(1) - nodes.at(0);
        Vec3 ac = nodes.at(2) - nodes.at(0);
        return 0.5 * ab.cross(ac);
    }

    Vec3 Triangle::calc_centroid() const
    {
        assert(nodes.size() == 3);
        return (nodes[0] + nodes[1] + nodes[2]) / 3; // 1/3(a + b + c)
    }

    Scalar Tetrahedron::calc_volume() const
    {
        assert(nodes.size() == 4);
        // Using vol = 1/3 * area_base * height = 1/3 * < 0.5 *  ab X ac , ad>
        Vec3 ab = nodes.at(1) - nodes.at(0);
        Vec3 ac = nodes.at(2) - nodes.at(0);
        Vec3 ad = nodes.at(3) - nodes.at(0);
        Scalar vol = -ab.cross(ac).dot(ad) / 6;
        assert(vol > 0);
        return vol;
    }

    Vec3 Tetrahedron::calc_centroid() const
    {
        assert(nodes.size() == 4);
        return 0.25 * (nodes[0] + nodes[1] + nodes[2] + nodes[3]); // 1/4(a + b + c + d)
    }

    void assign_face_properties(Vec3 &normal_area,
                                Vec3 &centroid_to_face_i,
                                Vec3 &centroid_to_face_j,
                                const Facegeom &face_geom,
                                const Vec3 &cell_center_i,
                                const Vec3 &cell_center_j)
    {
        Vec3 S_ij = face_geom.calc_area_normal();

        Scalar normal_dot_product = S_ij.dot(cell_center_j - cell_center_i);
        assert(normal_dot_product != 0); // Just banning this for now, altough it is possibly possible with a high skewness, but valid mesh
        if (normal_dot_product < 0)
            S_ij *= -1; // Flipping normal if it's not pointing from i to j
        normal_area = S_ij;
        Vec3 face_centroid = face_geom.calc_centroid();
        centroid_to_face_i = face_centroid - cell_center_i;
        centroid_to_face_j = face_centroid - cell_center_j;
    }

    Vec3 calc_ghost_centroid(Vec3 centroid_i, const Facegeom &boundary_face)
    {
        // distance from cell i to face V_i/ij = centroid_face - centroid_i
        // ghost centroid is located on opposite side of face with respect to i
        // centroid_ghost = centroid_i + 2*V_i/ij = 2*centroid_face - centroid_i
        Vec3 centroid_face = boundary_face.calc_centroid();
        return 2 * centroid_face - centroid_i;
    }
}

string Timer::get_elapsed_time() const
{
    using std::chrono::duration_cast;

    Time end_time = Clock::now();
    size_t total_miliseconds = duration_cast<Milliseconds>(end_time - start_time).count();
    size_t hours = total_miliseconds / (1000 * 60 * 60);
    size_t minutes = (total_miliseconds / (1000 * 60)) % 60;
    size_t seconds = (total_miliseconds / 1000) % 60;
    size_t milliseconds = total_miliseconds & 1000;

    using namespace std;
    constexpr size_t WIDTH = 4;
    std::stringstream ss;
    ss << "\n"
       << "    Hours:" << setw(WIDTH + 2) << hours << endl
       << "    Minutes:" << setw(WIDTH) << minutes << endl
       << "    Seconds:" << setw(WIDTH) << seconds << "." << milliseconds << endl;

    return ss.str();
};