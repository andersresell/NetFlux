#pragma once
#include "includes.hpp"


constexpr ShortIndex N_DIM{3};
constexpr ShortIndex N_VARS{N_DIM + 2};
constexpr ShortIndex N_TET_NODES{4};
constexpr ShortIndex N_TRI_NODES{3};
constexpr ShortIndex N_TET_FACES{4};

constexpr ShortIndex N_EQS_EULER{N_DIM + 2};

/*Generic container used for building other data types*/
template<typename T, Index N>
class Container{

protected:
    T data[N];
public:
    Container() : data{} {}

    Container(std::initializer_list<T> init){
        assert(init.size() == N);
        std::copy(init.begin(), init.end(), data);
    }
    
    T& operator[](Index i) {
        assert(i<N);
        return data[i];
    }
    const T& operator[](Index i) const{
        assert(i<N);
        return data[i];
    }
    
    Container<T,N> operator+(Container<T,N> rhs) const {
        for (Index i{0}; i<N; i++) rhs.data[i] += data[i];
        return rhs;
    };

    Container<T,N> operator-(const Container<T,N>& rhs) const {
        Container lhs{*this};
        for (Index i{0}; i<N; i++) lhs.data[i] -= rhs.data[i];
        return lhs;
    };

    void operator-=(const Container<T,N>& rhs) {
        for (Index i{0}; i<N; i++) data[i] -= rhs.data[i];
    };

    void operator+=(const Container<T,N>& rhs) {
        for (Index i{0}; i<N; i++) data[i] += rhs.data[i];
    };

    Container<T,N> operator*(T rhs) const {
        Container lhs{*this};
        for (Index i{0}; i<N; i++) lhs.data[i] *= rhs;
        return lhs;
    };

    void operator*=(T rhs){
        for (Index i{0}; i<N; i++) data[i] *= rhs;
    }

    Index size() const {return N;}
    
    friend std::ostream& operator<<(std::ostream& os, const Container& rhs){
        for (Index i{0}; i<N; i++) os << rhs[i] << " ";
        return os << std::endl;
    }
};


// struct ConsVars final : public Container<double, N_EQS_EULER>{
//     double rho() const {return data[0];}
//     double rho_u() const {return data[1];}
//     double rho_v() const {return data[2];}
//     double rho_w() const {return data[3];}
//     double E() const {return data[4];}
// };

// struct PrimVars final : public Container<double, N_EQS_EULER>{
//     double rho() const {return data[0];}
//     double u() const {return data[1];}
//     double v() const {return data[2];}
//     double w() const {return data[3];}
//     double p() const {return data[4];}
// };

namespace flow{
    
    /*template<Index N_EQS>
    struct FlowVar : public Container<double, N_EQS> {
        FlowVar() = default;
        FlowVar(std::initializer_list<double> init) : Container(init) {} 
    };*/

    constexpr double GAMMA{1.4};
    constexpr double GAMMA_MINUS_ONE{1-GAMMA};
    constexpr double GAMMA_MINUS_ONE_INV{1/GAMMA_MINUS_ONE};

    /*Variables of the Euler / NS equations, can be both conservative, primitive, fluxes, etc.*/
    struct EulerVar : public Container<double, N_EQS_EULER>{
        EulerVar() =  default;
        EulerVar(std::initializer_list<double> init) : Container(init) {}

        static EulerVar prim_to_cons(const EulerVar& V);
        static EulerVar cons_to_prim(const EulerVar& U);
        static double pressure(const EulerVar& U);

    };

    /*perhaps for later*/
    struct NS_VAR : public EulerVar{

    };
}
// --------------------------------------------------------------------
// Some enums that specify solver behaviour
// --------------------------------------------------------------------

enum class TimeScheme{ExplicitEuler, TVD_RK3};

enum class GradientScheme{Gauss};

enum class ConvScheme{Rusanov};

enum class Limiter{Minmod};

enum class BoundaryType{NoSlipWall, SlipWall, FarField};

const map<string, BoundaryType> boundary_type_from_string{
    {"NoSlipWall", BoundaryType::NoSlipWall},
    {"SlipWall", BoundaryType::SlipWall},
    {"FarField", BoundaryType::FarField},
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
        Vector<Index> boundary_face_indices;
    };

    /*Connectivity of a tetrahedron*/
    struct TetConnect final : public Container<Index, N_TET_NODES> {
        TetConnect() =  default;
        TetConnect(std::initializer_list<Index> init) : Container{init} {}
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
    struct TriConnect final : public Container<Index, N_TRI_NODES> {
        TriConnect() = default;
        TriConnect(std::initializer_list<Index> init) : Container{init} {}
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