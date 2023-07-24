#pragma once
#include "../containers/StaticContainer.hpp"
#include "../Config.hpp"

namespace geometry
{
    /*--------------------------------------------------------------------
    PrimalGrid stores reads the mesh file and stores nodes and
    connectivities of the FE type native/primal mesh, which is used to
    construct the face-based structure used by the FV solver.
    --------------------------------------------------------------------*/
    class PrimalGrid
    {
        friend class FV_Grid;
        /*These are read from the input mesh*/
        Vector<Vec3> nodes;
        Elements volume_elements;
        Vector<ElementPatch> element_patches;

        Elements face_elements; /*Storing elements of all faces*/

        void read_mesh(const Config &config);
        void read_netflux_mesh(const Config &config);
        void read_su2_mesh(const Config &config);

    public:
        const Vector<Vec3> &get_nodes() const { return nodes; }
        const Elements &get_volume_elements() const { return volume_elements; }
        const Elements &get_face_elements() const { return face_elements; }
        const Vector<ElementPatch> &get_element_patches() const { return element_patches; }
        Elements &get_face_elements() { return face_elements; }
    };

    /*Connectivities of elements laid out on the compressed row storage format (CRS) as METIS uses.
    It handles arbitrary elements.*/

    enum class ElementType
    {
        Triangle,
        Tetrahedron,
        Quadrilateral,
        Hexahedron,
        Pyramid
    };

    constexpr ShortIndex N_NODES_TET{4};
    constexpr ShortIndex N_FACES_TET{4};

    constexpr ShortIndex N_NODES_TRI{3};

    constexpr ShortIndex N_NODES_HEX{8};
    constexpr ShortIndex N_FACES_HEX{6};

    constexpr ShortIndex N_NODES_QUAD{4};

    constexpr ShortIndex N_NODES_PYRAMID{5};

    /*These constants store the maximum number of nodes among all element types
    for volume and face elements*/
    constexpr ShortIndex MAX_NODES_VOLUME_ELEMENT{N_NODES_HEX};
    constexpr ShortIndex MAX_NODES_FACE_ELEMENT{N_NODES_QUAD};

    const map<ElementType, ShortIndex> num_nodes_in_element = {{ElementType::Triangle, N_NODES_TRI},
                                                               {ElementType::Quadrilateral, N_NODES_QUAD},
                                                               {ElementType::Tetrahedron, N_NODES_TET},
                                                               {ElementType::Hexahedron, N_NODES_HEX},
                                                               {ElementType::Pyramid, N_NODES_PYRAMID}};
    const map<ElementType, string> element_type_to_string = {{ElementType::Triangle, "Triangle"},
                                                             {ElementType::Quadrilateral, "Quadrilateral"},
                                                             {ElementType::Tetrahedron, "Terahedron"},
                                                             {ElementType::Hexahedron, "Hexahedron"},
                                                             {ElementType::Pyramid, "Pyramid"}};
    const map<ElementType, bool> is_volume_element = {{ElementType::Triangle, false},
                                                      {ElementType::Tetrahedron, true},
                                                      {ElementType::Quadrilateral, false},
                                                      {ElementType::Hexahedron, true}};
    /*Identifiers used in the su2 mesh format: https://su2code.github.io/docs/Mesh-File/
    Add additional types when imlementing more elements*/
    enum class SU2_ElementID
    {
        Triangle = 5,
        Quadrilateral = 9,
        Tetrahedron = 10,
        Hexahedron = 12,
        Pyramid = 14
    };
    const map<SU2_ElementID, ElementType> su2_identifier_to_element_type = {{SU2_ElementID::Triangle, ElementType::Triangle},
                                                                            {SU2_ElementID::Quadrilateral, ElementType::Quadrilateral},
                                                                            {SU2_ElementID::Tetrahedron, ElementType::Tetrahedron},
                                                                            {SU2_ElementID::Hexahedron, ElementType::Hexahedron},
                                                                            {SU2_ElementID::Pyramid, ElementType::Pyramid}};
    class Elements
    {
    protected:
        Vector<Index> n_ptr = {0};
        Vector<Index> n_ind;
        Vector<ElementType> element_types;

    public:
        Elements() { assert(n_ptr.size() == 1); }

        Index size() const { return n_ptr.size() - 1; }

        void add_element(ElementType e_type, const Index *element)
        {
            ShortIndex n_nodes = num_nodes_in_element.at(e_type);
            n_ptr.emplace_back(n_ptr.back() + n_nodes);
            for (ShortIndex i{0}; i < n_nodes; i++)
                n_ind.emplace_back(element[i]);
            element_types.emplace_back(e_type);
        }

        const Index *get_element_nodes(Index i) const
        {
            assert((n_ptr[i + 1] - n_ptr[i]) == num_nodes_in_element.at(element_types[i]));
            return &n_ind[n_ptr[i]];
        }

        ElementType get_element_type(Index i) const
        {
            return element_types[i];
        }

        void reserve(Index n_elements, ShortIndex max_nodes_per_element);
        void shrink_to_fit();
    };

    // Vec3 calc_element_centroid(ElementType type, const Index *element, const Vector<Vec3> &nodes);

    // /*Connectivity of a tetrahedron*/
    // struct TetConnect final : public StaticContainer1D<Index, N_TET_NODES>
    // {
    //     TetConnect() = default;
    //     TetConnect(Index a, Index b, Index c, Index d) : StaticContainer1D{a, b, c, d} {}
    //     Index &a() { return data[0]; }
    //     Index &b() { return data[1]; }
    //     Index &c() { return data[2]; }
    //     Index &d() { return data[3]; }
    //     Index a() const { return data[0]; }
    //     Index b() const { return data[1]; }
    //     Index c() const { return data[2]; }
    //     Index d() const { return data[3]; }
    // };
    // /*Connectivity of a triangle*/
    // struct TriConnect : public StaticContainer1D<Index, N_TRI_NODES>
    // {
    //     TriConnect() = default;
    //     TriConnect(Index a, Index b, Index c) : StaticContainer1D{a, b, c} {}
    //     Index &a() { return data[0]; }
    //     Index &b() { return data[1]; }
    //     Index &c() { return data[2]; }
    //     Index a() const { return data[0]; }
    //     Index b() const { return data[1]; }
    //     Index c() const { return data[2]; }
    // };

    // struct SortedTriConnect : public TriConnect
    // {
    //     SortedTriConnect(const TriConnect &tc) : TriConnect{tc}
    //     {
    //         this->sort();
    //     }
    //     bool operator<(const SortedTriConnect &rhs) const
    //     {
    //         for (ShortIndex i{0}; i < N_TRI_NODES; i++)
    //             if (data[i] != rhs.data[i])
    //                 return data[i] < rhs.data[i];
    //         return false;
    //     }
    // };

    struct FaceElement
    {
        std::array<Index, MAX_NODES_FACE_ELEMENT> nodes;
        const ShortIndex n_nodes;
        const ElementType e_type;
        FaceElement(ElementType e_type, const Index *element) : n_nodes{num_nodes_in_element.at(e_type)}, e_type{e_type}
        {
            assert(is_volume_element.at(e_type));
            std::copy(element, element + n_nodes, nodes);
        }
    };

    inline FaceElement get_face_element_k_of_volume_element(ElementType volume_element_type,
                                                            const Index *volume_element,
                                                            ShortIndex face_k);
    inline void get_face_element_k_of_tetrahedron(ElementType volume_element_type,
                                                  const Index *ve,
                                                  ShortIndex face_k,
                                                  ElementType &face_element_type,
                                                  array<Index, MAX_NODES_FACE_ELEMENT> &fe);
    inline void get_face_element_k_of_hexahedron(ElementType volume_element_type,
                                                 const Index *ve,
                                                 ShortIndex face_k,
                                                 ElementType &face_element_type,
                                                 array<Index, MAX_NODES_FACE_ELEMENT> &fe);

    struct SortedFaceElement : public FaceElement
    {
        SortedFaceElement(const FaceElement &rhs) : FaceElement(rhs)
        {
            std::sort(nodes.begin(), nodes.begin() + n_nodes);
        }

        bool operator<(const SortedFaceElement &other) const
        {
            if (n_nodes == other.n_nodes)
            {
                for (ShortIndex i{0}; i < n_nodes; i++)
                    if (nodes[i] != other.nodes[i])
                        return nodes[i] < other.nodes[i];
                return false;
            }
            else
                return n_nodes < other.n_nodes;
        }
    };

    /*Holds name and triangles of a boundary patch*/
    struct ElementPatch
    {
        string patch_name;
        Elements boundary_elements;
    };

    // returns the node connectivity of face_k from the node connectivity of tetraheder
    // TriConnect tet_face_connectivity(TetConnect tc, ShortIndex face_k);

    /*--------------------------------------------------------------------
    Function to calculate geometry properties of various elements
    --------------------------------------------------------------------*/
    inline void volume_element_calc_geometry_properties(ElementType e_type,
                                                        const Index *element,
                                                        const Vector<Vec3> &nodes,
                                                        Scalar &volume,
                                                        Vec3 &centroid);

    inline void face_element_calc_centroid(ElementType e_type,
                                           const Index *element,
                                           const Vector<Vec3> &nodes,
                                           Vec3 &centroid);
    inline void face_element_calc_face_normal(ElementType e_type, const Index *element, const Vector<Vec3> &nodes, Vec3 &S_ij);
    inline void tetrahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                                     Scalar &volume, Vec3 &centroid);
    inline void tetrahedron_calc_volume(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Scalar &volume);
    inline void tetrahedron_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid);

    inline void triangle_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &S_ij);
    inline void triangle_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &centroid);

    inline void pyramid_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, const Vec3 &n4,
                                                 Scalar &volume, Vec3 &centroid);

    inline void hexahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                                    const Vec3 &n4, const Vec3 &n5, const Vec3 &n6, const Vec3 &n7,
                                                    Scalar &volume, Vec3 &centroid);
    inline void quadrilateral_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &S_ij);
    inline void quadrilateral_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid);

    // /*Astract face geometry class*/
    // struct Facegeom
    // {
    //     Vector<Vec3> nodes;
    //     Facegeom(){};
    //     virtual Vec3 calc_area_normal() const = 0;
    //     virtual Vec3 calc_centroid() const = 0;
    // };

    // struct Triangle final : Facegeom
    // {
    //     Triangle(Vec3 a, Vec3 b, Vec3 c) : Facegeom() { nodes = {a, b, c}; }
    //     Vec3 calc_area_normal() const final;
    //     Vec3 calc_centroid() const final;
    // };

    // /*Abstract volume geometry class*/
    // struct Polyhedra
    // {
    //     Vector<Vec3> nodes;
    //     virtual Scalar calc_volume() const = 0;
    //     virtual Vec3 calc_centroid() const = 0;
    // };

    // struct Tetrahedron final : Polyhedra
    // {
    //     Tetrahedron(Vec3 a, Vec3 b, Vec3 c, Vec3 d) : Polyhedra() { nodes = {a, b, c, d}; }
    //     Scalar calc_volume() const final;
    //     Vec3 calc_centroid() const final;
    // };

}