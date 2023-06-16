#pragma once
#include "includes.hpp"
#include "containers/StaticContainer.hpp"

constexpr ShortIndex N_TET_NODES{4};
constexpr ShortIndex N_TRI_NODES{3};
constexpr ShortIndex N_TET_FACES{4};

constexpr ShortIndex N_EQS_EULER{N_DIM + 2};



// --------------------------------------------------------------------
// Some enums that specify solver behaviour and corresponding string -> enum maps for some
// --------------------------------------------------------------------
enum class MainSolverType{Euler}; 

    const map<string, MainSolverType> main_solver_from_string{
        {"Euler", MainSolverType::Euler}
    };

enum class SolverType{Euler}; //The main solver can (but need not) be comprised of multiple sub solvers

enum class TimeIntegrationType{Explicit, Implicit};

enum class TimeScheme{ExplicitEuler, TVD_RK3};

    const map<string, TimeScheme> time_scheme_from_string{
        {"ExplicitEuler", TimeScheme::ExplicitEuler},
        {"TVD_RK3", TimeScheme::TVD_RK3}
    };


enum class SpatialOrder{First, Second};

    const map<string, SpatialOrder> spatial_order_from_string{
        {"First", SpatialOrder::First},
        {"Second", SpatialOrder::Second}
    };

enum class GradientScheme{GreenGauss};

    const map<string, GradientScheme> gradient_scheme_from_string{
        {"GreenGauss", GradientScheme::GreenGauss}
    };

enum class InviscidFluxScheme{Rusanov, HLLC};

    const map<string, InviscidFluxScheme>  inviscid_flux_scheme_from_string{
        {"Rusanov", InviscidFluxScheme::Rusanov},
        {"HLLC", InviscidFluxScheme::HLLC}
    };

enum class Limiter{NONE, Barth};

    const map<string, Limiter>  limiter_from_string{
        {"NONE", Limiter::NONE},
        {"Barth", Limiter::Barth}
    };

enum class BoundaryType{NoSlipWall, SlipWall, FarField};

    const map<string, BoundaryType>  boundary_type_from_string{
        {"NoSlipWall", BoundaryType::NoSlipWall},
        {"SlipWall", BoundaryType::SlipWall},
        {"FarField", BoundaryType::FarField}
    };

enum class InitialConditionOption{Freestream};

    const map<string, InitialConditionOption>  initial_condition_option_from_string{
        {"Freestream", InitialConditionOption::Freestream}
    };




namespace geom{

    struct Face{
        Face(Index i, Index j, Vec3 S_ij, Vec3 r_im, Vec3 r_jm) : S_ij{S_ij}, i{i}, j{j}, r_im{r_im}, r_jm{r_jm} {}
        Face(Index i, Index j) : i{i}, j{j} {} 
        Vec3 S_ij; //Area normal vector from cell i to j
        Index i,j; //Indices of cell i and j
        Vec3 r_im, r_jm; //vectors from each cell center to the face centroid
        friend std::ostream& operator<<(std::ostream& os, const Face& f) {
            os << "(i,j) = ("<< f.i << "," << f.j << "), S_ij: " + horizontal_string_Vec3(f.S_ij) +
            ", r_im: " + horizontal_string_Vec3(f.r_im) + ", r_jm: " + horizontal_string_Vec3(f.r_jm) << endl;
            return os;
        }

        bool operator<(const Face& rhs) const {
            if (i != rhs.i) 
                return i < rhs.i;
    
            assert(j != rhs.j); //This means cell indices are identical 
                
            return j < rhs.j;
        }

    };

    struct Cell{
        Cell() = default;
        Cell(double cell_volume, Vec3 centroid) : cell_volume{cell_volume}, centroid{centroid} {}
        double cell_volume;
        Vec3 centroid;
        friend std::ostream& operator<<(std::ostream& os, const Cell& c) {
            os << "vol = " << c.cell_volume << ", centroid = " + horizontal_string_Vec3(c.centroid) << endl; 
            return os;
        }
    };


    struct Patch{
        BoundaryType boundary_type;
        //Vector<Index> boundary_face_indices;
        Index N_FACES;
        Index FIRST_FACE;
    };

    /*Connectivity of a tetrahedron*/
    struct TetConnect final : public StaticContainer1D<Index, N_TET_NODES> {
        TetConnect() =  default;
        TetConnect(Index a, Index b, Index c, Index d) : StaticContainer1D{a,b,c,d} {}
        Index& a() {return data[0];}
        Index& b() {return data[1];}
        Index& c() {return data[2];}
        Index& d() {return data[3];}
        Index a() const {return data[0];}
        Index b() const {return data[1];}
        Index c() const {return data[2];}
        Index d() const {return data[3];}
    };
    /*Connectivity of a triangle*/
    struct TriConnect final : public StaticContainer1D<Index, N_TRI_NODES> {
        TriConnect() = default;
        TriConnect(std::initializer_list<Index> init) : StaticContainer1D{init} {}
        Index& a() {return data[0];}
        Index& b() {return data[1];}
        Index& c() {return data[2];}
        Index a() const {return data[0];}
        Index b() const {return data[1];}
        Index c() const {return data[2];}
    };
    /*Holds name and triangles of a boundary patch*/
    struct TriPatchConnect{
        string patch_name;
        Vector<TriConnect> triangles;
    };

    //returns the node connectivity of face_i from the node connectivity of tetraheder  
    TriConnect tet_face_connectivity(TetConnect tc, ShortIndex face_k); 

    /*Astract face geometry class*/
    struct Facegeom{
        Vector<Vec3> nodes;
        Facegeom() {};
        virtual Vec3 calc_area_normal() const = 0;
        virtual Vec3 calc_centroid() const = 0;
    };

    struct Triangle final : Facegeom{
        Triangle(Vec3 a, Vec3 b, Vec3 c) : Facegeom() { nodes = {a,b,c};}
        Vec3 calc_area_normal() const final;
        Vec3 calc_centroid() const final;
    };

    /*Abstract volume geometry class*/
    struct Polyhedra{
        Vector<Vec3> nodes;
        virtual double calc_volume() const = 0;
        virtual Vec3 calc_centroid() const = 0;
    };

    struct Tetrahedron final : Polyhedra{
        Tetrahedron(Vec3 a, Vec3 b, Vec3 c, Vec3 d) : Polyhedra() {nodes = {a,b,c,d};}
        double calc_volume() const final;
        Vec3 calc_centroid() const final;

    };

    void assign_face_properties(Face& face, const Facegeom& face_geom, const Vec3& cell_center_i, const Vec3& cell_center_j);

    Vec3 calc_ghost_centroid(Vec3 centroid_i, const Facegeom& boundary_face);
}