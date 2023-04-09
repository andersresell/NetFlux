#include "../include/Utilities.hpp"


namespace Geometry{
        
    
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
}
