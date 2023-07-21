#pragma once
#include "../containers/StaticContainer.hpp"
#include "../Config.hpp"

namespace geometry
{
    struct Connectivities;
    /*--------------------------------------------------------------------
    PrimalGrid stores reads the mesh file and stores nodes and
    connectivities of the FE type native/primal mesh, which is used to
    construct the face-based structure used by the FV solver.
    --------------------------------------------------------------------*/
    class PrimalGrid
    {

        Vector<Vec3> nodes;
        Elements elements;
        Vector<PatchElements> patches;

        void read_mesh(const Config &config);
        void read_netflux_mesh(const Config &config);
        void read_su2_mesh(const Config &config);

    public:
        const Vector<Vec3> &get_nodes() const { return nodes; }
        const Elements &get_elements() const { return elements; }
    };

    /*Connectivities of elements laid out on the compressed row storage format (CRS) as METIS uses.
    It handles arbitrary elements.*/

    enum class ElementType
    {
        Tri,
        Tet,
        Quad,
        Hex
    };

    const map<ElementType, ShortIndex> num_nodes_in_element = {{ElementType::Tri, 3},
                                                               {ElementType::Quad, 4},
                                                               {ElementType::Tet, 4},
                                                               {ElementType::Hex, 8}};
    const map<ElementType, bool> is_volume_element = {{ElementType::Tri, false},
                                                      {ElementType::Tet, true},
                                                      {ElementType::Quad, false},
                                                      {ElementType::Hex, true}};
    class Elements
    {
        Vector<Index> n_ptr;
        Vector<Index> n_ind;
        Vector<ElementType> element_types;

    public:
        Index size() const { return n_ptr.size() - 1; }

        const Index *get_element_nodes(Index i) const
        {
            assert((n_ptr[i + 1] - n_ptr[i]) == num_nodes_in_element.at(element_types[i]));
            return &n_ind[n_ptr[i]];
        }

        ElementType get_element_type(Index i) const
        {
            return element_types[i];
        }
    };

    // Vec3 calc_element_centroid(ElementType type, const Index *element, const Vector<Vec3> &nodes);

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
    struct PatchElements
    {
        string patch_name;
        Elements surface_elements;
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