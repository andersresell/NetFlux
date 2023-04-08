#pragma once
#include "includes.hpp"
#include <stdint.h>


using Index = uint32_t; //Used for accessing cell indices, etc 
using ShortIndex = uint8_t; //Used for looping over dimensions, etc

constexpr ShortIndex N_DIM{3};
constexpr ShortIndex N_VARS{N_DIM + 2};

using Vec3 = Eigen::Vector3d;

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
    Face(Vec3 S_ij, Index i, Index j) : S_ij{S_ij}, i{i}, j{j} {}
    Vec3 S_ij; //Area normal vector from cell i to j
    Index i,j; //Indices of cell i and j
};

struct Cell{
    double cell_volume;
    Vec3 centroid;
};

using FaceContainer = vector<Face>;
using CellContainer = vector<Cell>;


enum class BoundaryType{NoSlipWall, SlipWall, FarField};

const map<string, BoundaryType> boundary_type_from_string{
    {"NoSlipWall", BoundaryType::NoSlipWall},
    {"SlipWall", BoundaryType::SlipWall},
    {"FarField", BoundaryType::FarField},
};
