#pragma once
#include "includes.hpp"
#include <stdint.h>


using Index = uint32_t; //Used for accessing cell indices, etc 
using ShortIndex = uint8_t; //Used for looping over dimensions, etc

constexpr ShortIndex N_DIM{3};
constexpr ShortIndex N_VARS{N_DIM + 2};
constexpr ShortIndex N_TET_NODES{4};
constexpr ShortIndex N_TRI_NODES{3};
constexpr ShortIndex N_TET_FACES{4};

using Vec3 = Eigen::Vector3d;

inline string horizontal_string_Vec3(const Vec3& v) {return "["+std::to_string(v[0])+", "+std::to_string(v[1])+", "+std::to_string(v[2])+" ]";}

template<typename T, size_t N>
inline bool arrays_equal(const array<T, N>& a,const array<T, N>& b){
    assert (a.size() == b.size());
    for (size_t i{0}; i<a.size(); i++) 
        if (a[i] != b[i]) return false;
    return true;
}

struct FlowVec5{
    double u1, u2, u3, u4, u5;

        FlowVec5 operator+(const FlowVec5 &rhs) const { return {u1 + rhs.u1, u2 + rhs.u2, u3 + rhs.u3, u4 + rhs.u4, u5 + rhs.u5}; }

        FlowVec5 operator-(const FlowVec5 &rhs) const { return {u1 - rhs.u1, u2 - rhs.u2, u3 - rhs.u3, u4 - rhs.u4, u5 + rhs.u5}; }

        void operator*=(double rhs) {
            u1 *= rhs;
            u2 *= rhs;
            u3 *= rhs;
            u4 *= rhs;
            u5 *= rhs;
        }
};

struct ConservativeVars final : public FlowVec5{
    double rho() const {return u1;}
    double rho_u() const {return u2;}
    double rho_v() const {return u3;}
    double rho_w() const {return u4;}
    double E() const {return u5;}
  
};

struct PrimitiveVars final : public FlowVec5{
    double rho() const {return u1;}
    double u() const {return u2;}
    double v() const {return u3;}
    double w() const {return u4;}
    double p() const {return u5;}
};


using SolutionContainer = vector<ConservativeVars>;
using FluxContainer = vector<FlowVec5>;

struct Face{
    Face(Index i, Index j, Vec3 S_ij, Vec3 r_im, Vec3 r_jm) : S_ij{S_ij}, i{i}, j{j}, r_im{r_im}, r_jm{r_jm} {}
    Face(Index i, Index j) : i{i}, j{j} {} 
    Vec3 S_ij; //Area normal vector from cell i to j
    Index i,j; //Indices of cell i and j
    Vec3 r_im, r_jm; //vectors from each cell center to the face centroid
    friend std::ostream& operator<<(std::ostream& os, const Face& f) {
        os << "i: "<< f.i << ", j: " << f.j << ", S_ij: " + horizontal_string_Vec3(f.S_ij) +
        ", r_im: " + horizontal_string_Vec3(f.r_im) + ", r_jm: " + horizontal_string_Vec3(f.r_jm) << endl;
        return os;
    }
};

struct Cell{
    Cell() {}
    Cell(double cell_volume, Vec3 centroid) : cell_volume{cell_volume}, centroid{centroid} {}
    double cell_volume;
    Vec3 centroid;
    friend std::ostream& operator<<(std::ostream& os, const Cell& c) {
        os << "cell vol = " << c.cell_volume << ", centroid = " + horizontal_string_Vec3(c.centroid) << endl; 
        return os;
    }
};

using FaceContainer = vector<Face>;
using CellContainer = vector<Cell>;


enum class BoundaryType{NoSlipWall, SlipWall, FarField};

const map<string, BoundaryType> boundary_type_from_string{
    {"NoSlipWall", BoundaryType::NoSlipWall},
    {"SlipWall", BoundaryType::SlipWall},
    {"FarField", BoundaryType::FarField},
};

struct Patch{
    BoundaryType boundary_type;
    vector<Index> boundary_face_indices;
};

namespace Geometry{
    struct TetConnectivity{
        Index a,b,c,d;
    };
    struct TriConnectivity{
        Index a,b,c;
    };
    struct TriPatchConnectivity{
        BoundaryType boundary_type;
        vector<TriConnectivity> triangles;
    };

    //returns the node connectivity of face_i from the node connectivity of tetraheder  
    TriConnectivity tet_face_connectivity(TetConnectivity tc, ShortIndex face_k); 

    struct FaceGeom{
        vector<Vec3> nodes;
        FaceGeom() {};
        virtual Vec3 calc_area_normal() const = 0;
        virtual Vec3 calc_centroid() const = 0;
    };

    struct Triangle final : FaceGeom{
        Triangle(Vec3 a, Vec3 b, Vec3 c) : FaceGeom() { nodes = {a,b,c};}
        Vec3 calc_area_normal() const final;
        Vec3 calc_centroid() const final;
    };

    struct Polyhedra{
        vector<Vec3> nodes;
        virtual double calc_volume() const = 0;
        virtual Vec3 calc_centroid() const = 0;
    };

    struct Tetrahedron final : Polyhedra{
        Tetrahedron(Vec3 a, Vec3 b, Vec3 c, Vec3 d) : Polyhedra() {nodes = {a,b,c,d};}
        double calc_volume() const final;
        Vec3 calc_centroid() const final;

    };

    void assign_face_properties(Face& face, const FaceGeom& face_geom, const Vec3& cell_center_i, const Vec3& cell_center_j);

    
}