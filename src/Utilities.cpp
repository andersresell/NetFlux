#include "../include/Utilities.hpp"

namespace flow{
    
    FlowVar FlowVar::prim_to_cons(const FlowVar& V){
        return {V[0], 
                V[0]*V[1],
                V[0]*V[2],
                V[0]*V[3],
                GAMMA_MINUS_ONE_INV * V[4] + 0.5*V[0]*(V[1]*V[1] + V[2]*V[2] + V[3]*V[3])};
    }
    FlowVar FlowVar::cons_to_prim(const FlowVar& U){
        return {U[0], 
                U[1]/U[0],
                U[2]/U[0],
                U[3]/U[0],
                pressure(U)};
        }    

    
    double FlowVar::pressure(const FlowVar& U){
        return GAMMA_MINUS_ONE * (U[4] - 0.5/U[0]*(U[1]*U[1] + U[2]*U[2] + U[3]*U[3]));
    }
}




namespace Geom{
        
    TriConnect tet_face_connectivity(TetConnect tc, ShortIndex face_k) {
        assert(face_k >= 0 && face_k < N_TET_FACES);
        switch (face_k){
            case 0:
                return {tc.a(), tc.b(), tc.c()};
            case 1:
                return {tc.a(), tc.b(), tc.d()};
            case 2:
                return {tc.a(), tc.c(), tc.d()};
            case 3:
                return {tc.b(), tc.c(), tc.d()};
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
        double vol = -ab.cross(ac).dot(ad) / 6; 
        assert(vol > 0);
        return vol;
    }

    Vec3 Tetrahedron::calc_centroid() const {
        assert(nodes.size() == 4);
        return 0.25 * (nodes[0] + nodes[1] + nodes[2] + nodes[3]); // 1/4(a + b + c + d)
    }

    Face create_face_from_geom(Index i, Index j, const FaceGeom& face_geom, Vec3 cell_center_i, Vec3 cell_center_j){
        Vec3 face_centroid = face_geom.calc_centroid();
        return {i, j, face_geom.calc_area_normal(), face_centroid - cell_center_i, face_centroid - cell_center_j};
    }

    void assign_face_properties(Face& face, const FaceGeom& face_geom, const Vec3& cell_center_i, const Vec3& cell_center_j){
        Vec3 S_ij = face_geom.calc_area_normal();
        
        double normal_dot_product = S_ij.dot(cell_center_j - cell_center_i); 
        assert(normal_dot_product != 0); //Just banning this for now, altough it is possibly possible with a high skewness, but valid mesh
        if (normal_dot_product < 0) S_ij *= -1; //Flipping normal if it's not pointing from i to j
        face.S_ij = S_ij;
        Vec3 face_centroid = face_geom.calc_centroid();
        face.r_im = cell_center_i - face_centroid;  
        face.r_jm = cell_center_j - face_centroid;
    }
    Vec3 calc_ghost_centroid(Vec3 centroid_i, const FaceGeom& boundary_face){
        //distance from cell i to face V_i/ij = centroid_face - centroid_i
        //ghost centroid is located on opposite side of face with respect to i
        //centroid_ghost = centroid_i + 2*V_i/ij = 2*centroid_face - centroid_i
        Vec3 centroid_face = boundary_face.calc_centroid();
        return 2 * centroid_face - centroid_i; 

        
    }
}
