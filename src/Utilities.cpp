#include "../include/Utilities.hpp"


template<typename T, size_t N> bool arrays_equal(const array<T, N>& a,const array<T, N>& b){
    assert (a.size() == b.size());
    for (size_t i{0}; i<a.size(); i++) 
        if (a[i] != b[i]) return false;
    return true;
}

namespace Geometry{
        
    TriConnectivity tet_face_connectivity(TetConnectivity tc, ShortIndex face_k) {
        assert(face_k >= 0 && face_k < N_TET_FACES);
        switch (face_k){
            case 0:
                return {tc.a, tc.b, tc.c};
            case 1:
                return {tc.a, tc.b, tc.c};
            case 2:
                return {tc.a, tc.c, tc.d};
            case 3:
                return {tc.b, tc.c, tc.d};
        }
    }
    
    Vec3 Triangle::calc_area_normal() const {
        assert(nodes.size() == 3);
        //Using that the cross product of two vectors is the area of a parallelogram (with correct normal), so the half is a the area of a triangle
        Vec3 ab = nodes.at(1) - nodes.at(0);
        Vec3 ac = nodes.at(2) - nodes.at(0);
        return 0.5 * ab.cross(ac);
    }

    Vec3 Triangle::calc_centroid() const {
        assert(nodes.size() == 3);
        return (nodes[0] + nodes[1] + nodes[2])/3; // 1/3(a + b + c)
    }

    double Tetrahedron::calc_volume() const {
        assert(nodes.size() == 4);
        //Using vol = 1/3 * area_base * height = 1/3 * < 0.5 *  ab X ac , ad>
        Vec3 ab = nodes.at(1) - nodes.at(0);
        Vec3 ac = nodes.at(2) - nodes.at(0);
        Vec3 ad = nodes.at(3) - nodes.at(0);
        return ab.cross(ac).dot(ad) / 6;
    }

    Vec3 Tetrahedron::calc_centroid() const {
        assert(nodes.size() == 4);
        return 0.25 * (nodes[0] + nodes[1] + nodes[2] + nodes[3]); // 1/4(a + b + c + d)
    }


    void assign_face_properties(Face& face, const FaceGeom& face_geom, const Vec3& cell_center_i, const Vec3& cell_center_j){
        face.S_ij = face_geom.calc_area_normal();
        Vec3 face_centroid = face_geom.calc_centroid();
        face.r_im = cell_center_i - face_centroid;
        face.r_jm = cell_center_j - face_centroid;
    }

    Face create_face_from_geom(Index i, Index j, const FaceGeom& face_geom, Vec3 cell_center_i, Vec3 cell_center_j){
        Vec3 face_centroid = face_geom.calc_centroid();
        return {i, j, face_geom.calc_area_normal(), face_centroid - cell_center_i, face_centroid - cell_center_j};
    }


}
