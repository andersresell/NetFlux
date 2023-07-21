#pragma once
#include "Includes.hpp"
#include "containers/StaticContainer.hpp"

constexpr ShortIndex N_TET_NODES{4};
constexpr ShortIndex N_TRI_NODES{3};
constexpr ShortIndex N_TET_FACES{4};

constexpr ShortIndex N_EQS_EULER{N_DIM + 2};

// --------------------------------------------------------------------
// Some enums that specify solver behaviour and corresponding string -> enum maps for some
// --------------------------------------------------------------------
enum class MainSolverType
{
    Euler
};

const map<string, MainSolverType> main_solver_from_string{
    {"Euler", MainSolverType::Euler}};

enum class SolverType
{
    Euler
}; // The main solver can (but need not) be comprised of multiple sub solvers

const map<SolverType, string> string_from_solver_type{
    {SolverType::Euler, "Euler"}};

enum class TimeIntegrationType
{
    Explicit,
    Implicit
};

enum class TimeScheme
{
    ExplicitEuler,
    TVD_RK3
};

const map<string, TimeScheme> time_scheme_from_string{
    {"ExplicitEuler", TimeScheme::ExplicitEuler},
    {"TVD_RK3", TimeScheme::TVD_RK3}};

enum class SpatialOrder
{
    First,
    Second
};

const map<string, SpatialOrder> spatial_order_from_string{
    {"First", SpatialOrder::First},
    {"Second", SpatialOrder::Second}};

enum class GradientScheme
{
    GreenGauss
};

const map<string, GradientScheme> gradient_scheme_from_string{
    {"GreenGauss", GradientScheme::GreenGauss}};

enum class InviscidFluxScheme
{
    Rusanov,
    HLLC
};

const map<string, InviscidFluxScheme> inviscid_flux_scheme_from_string{
    {"Rusanov", InviscidFluxScheme::Rusanov},
    {"HLLC", InviscidFluxScheme::HLLC}};

enum class Limiter
{
    NONE,
    Barth
};

const map<string, Limiter> limiter_from_string{
    {"NONE", Limiter::NONE},
    {"Barth", Limiter::Barth}};

enum class BoundaryType
{
    NoSlipWall,
    SlipWall,
    FarField
};

const map<string, BoundaryType> boundary_type_from_string{
    {"NoSlipWall", BoundaryType::NoSlipWall},
    {"SlipWall", BoundaryType::SlipWall},
    {"FarField", BoundaryType::FarField}};

enum class InitialConditionOption
{
    Freestream
};

const map<string, InitialConditionOption> initial_condition_option_from_string{
    {"Freestream", InitialConditionOption::Freestream}};

namespace standard_air
{
    constexpr Scalar gas_constant{287.058};
    constexpr Scalar gamma{1.4};
    constexpr Scalar density{1.225};
    constexpr Scalar pressure{101325.0};
    constexpr Scalar temperature = pressure / (density * gas_constant);
}

namespace geom
{
    class Grid;

    /*A structure of arrays (SoA) containing the faces and their required properties*/
    class Faces
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

    // struct Cell
    // {
    //     Cell() = default;
    //     Cell(Scalar cell_volume, Vec3 centroid) : cell_volume{cell_volume}, centroid{centroid} {}
    //     Scalar cell_volume;
    //     Vec3 centroid;
    //     friend std::ostream &operator<<(std::ostream &os, const Cell &c)
    //     {
    //         os << "vol = " << c.cell_volume << ", centroid = " + horizontal_string_Vec3(c.centroid) << endl;
    //         return os;
    //     }
    // };

    struct Patch
    {
        BoundaryType boundary_type;
        // Vector<Index> boundary_face_indices;
        Index N_FACES;
        Index FIRST_FACE;
    };

    /*Connectivity of a tetrahedron*/
    struct TetConnect final : public StaticContainer1D<Index, N_TET_NODES>
    {
        TetConnect() = default;
        TetConnect(Index a, Index b, Index c, Index d) : StaticContainer1D{a, b, c, d} {}
        Index &a() { return data[0]; }
        Index &b() { return data[1]; }
        Index &c() { return data[2]; }
        Index &d() { return data[3]; }
        Index a() const { return data[0]; }
        Index b() const { return data[1]; }
        Index c() const { return data[2]; }
        Index d() const { return data[3]; }
    };
    /*Connectivity of a triangle*/
    struct TriConnect : public StaticContainer1D<Index, N_TRI_NODES>
    {
        TriConnect() = default;
        TriConnect(Index a, Index b, Index c) : StaticContainer1D{a, b, c} {}
        Index &a() { return data[0]; }
        Index &b() { return data[1]; }
        Index &c() { return data[2]; }
        Index a() const { return data[0]; }
        Index b() const { return data[1]; }
        Index c() const { return data[2]; }
    };

    struct SortedTriConnect : public TriConnect
    {
        SortedTriConnect(const TriConnect &tc) : TriConnect{tc}
        {
            this->sort();
        }
        bool operator<(const SortedTriConnect &rhs) const
        {
            for (ShortIndex i{0}; i < N_TRI_NODES; i++)
                if (data[i] != rhs.data[i])
                    return data[i] < rhs.data[i];
            return false;
        }
    };

    /*Holds name and triangles of a boundary patch*/
    struct TriPatchConnect
    {
        string patch_name;
        Vector<TriConnect> triangles;
    };

    // returns the node connectivity of face_k from the node connectivity of tetraheder
    TriConnect tet_face_connectivity(TetConnect tc, ShortIndex face_k);

    /*Astract face geometry class*/
    struct Facegeom
    {
        Vector<Vec3> nodes;
        Facegeom(){};
        virtual Vec3 calc_area_normal() const = 0;
        virtual Vec3 calc_centroid() const = 0;
    };

    struct Triangle final : Facegeom
    {
        Triangle(Vec3 a, Vec3 b, Vec3 c) : Facegeom() { nodes = {a, b, c}; }
        Vec3 calc_area_normal() const final;
        Vec3 calc_centroid() const final;
    };

    /*Abstract volume geometry class*/
    struct Polyhedra
    {
        Vector<Vec3> nodes;
        virtual Scalar calc_volume() const = 0;
        virtual Vec3 calc_centroid() const = 0;
    };

    struct Tetrahedron final : Polyhedra
    {
        Tetrahedron(Vec3 a, Vec3 b, Vec3 c, Vec3 d) : Polyhedra() { nodes = {a, b, c, d}; }
        Scalar calc_volume() const final;
        Vec3 calc_centroid() const final;
    };

    void assign_face_properties(Vec3 &normal_area,
                                Vec3 &centroid_to_face_i,
                                Vec3 &centroid_to_face_j,
                                const Facegeom &face_geom,
                                const Vec3 &cell_center_i,
                                const Vec3 &cell_center_j);

    Vec3 calc_ghost_centroid(Vec3 centroid_i, const Facegeom &boundary_face);
}

class Timer
{
    using Clock = std::chrono::high_resolution_clock;
    using Time = std::chrono::_V2::system_clock::time_point;
    using Milliseconds = std::chrono::milliseconds;

    Time start_time; // Used for measuring time

public:
    void start_counter() { start_time = Clock::now(); }
    string get_elapsed_time() const;
};