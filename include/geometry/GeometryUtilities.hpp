#pragma once
#include "../Utilities.hpp"

namespace geometry
{

    /*The different element types implemented (Uses vtk numbering system)*/
    enum class ElementType
    {
        Triangle = 5,
        Quadrilateral = 9,
        Tetrahedron = 10,
        Hexahedron = 12,
        Pyramid = 14
    };
    inline bool legal_vtk_element_identifier(ShortIndex vtk_id)
    {
        auto e_type = static_cast<ElementType>(vtk_id);
        switch (e_type)
        {
        case ElementType::Triangle:
        case ElementType::Quadrilateral:
        case ElementType::Tetrahedron:
        case ElementType::Hexahedron:
        case ElementType::Pyramid:
            return true;
        default:
            return false;
        }
    }

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
                                                             {ElementType::Tetrahedron, "Tetrahedron"},
                                                             {ElementType::Hexahedron, "Hexahedron"},
                                                             {ElementType::Pyramid, "Pyramid"}};
    const map<ElementType, bool> is_volume_element = {{ElementType::Triangle, false},
                                                      {ElementType::Tetrahedron, true},
                                                      {ElementType::Quadrilateral, false},
                                                      {ElementType::Hexahedron, true}};

    /*Connectivities of elements laid out on the compressed row storage format (CRS) as METIS uses.*/
    class Elements
    {
        friend class Faces;

    protected:
        Vector<Index> n_ptr = {0};
        Vector<Index> n_ind;
        Vector<ElementType> element_types;

    public:
        Elements() = default;
        Elements(const Elements &other);
        Elements &operator=(Elements other);
        Index size() const { return n_ptr.size() - 1; }

        // const Vector<Index> &get_n_ptr() const { return n_ptr; }
        // const Vector<Index> &get_n_ind() const { return n_ind; }

        void add_element(ElementType e_type, const Index *element)
        {
            ShortIndex n_nodes = num_nodes_in_element.at(e_type);
            n_ptr.emplace_back(n_ptr.back() + n_nodes);
            for (ShortIndex k{0}; k < n_nodes; k++)
                n_ind.emplace_back(element[k]);
            element_types.emplace_back(e_type);
        }

        Index sum_nodes_over_all_elements() const { return n_ptr.back(); }

        const Index *get_element_nodes(Index i) const
        {
            assert((n_ptr[i + 1] - n_ptr[i]) == num_nodes_in_element.at(element_types[i]));
            return &n_ind[n_ptr[i]];
        }
        const ShortIndex get_num_nodes_of_element(Index i) const
        {
            Index num_nodes{n_ptr[i + 1] - n_ptr[i]};
            assert(num_nodes == num_nodes_in_element.at(element_types[i]));
            return num_nodes;
        }

        ElementType get_element_type(Index i) const
        {
            return element_types[i];
        }

        void reserve(Index n_elements, ShortIndex max_nodes_per_element);
        void shrink_to_fit();

        string to_string(Index i) const;
    };

    struct FaceElement
    {
        std::array<Index, MAX_NODES_FACE_ELEMENT> nodes;
        const ShortIndex n_nodes;
        const ElementType e_type;
        FaceElement(ElementType e_type, const Index *element) : n_nodes{num_nodes_in_element.at(e_type)}, e_type{e_type}
        {
            assert(!is_volume_element.at(e_type));
            std::copy(element, element + n_nodes, nodes.begin());
        }
    };

    FaceElement get_face_element_k_of_volume_element(ElementType volume_element_type,
                                                     const Index *volume_element,
                                                     ShortIndex face_k);
    void get_face_element_k_of_tetrahedron(ElementType volume_element_type,
                                           const Index *ve,
                                           ShortIndex face_k,
                                           ElementType &face_element_type,
                                           array<Index, MAX_NODES_FACE_ELEMENT> &fe);
    void get_face_element_k_of_hexahedron(ElementType volume_element_type,
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

    /*--------------------------------------------------------------------
    Functions to calculate geometry properties of various elements
    --------------------------------------------------------------------*/
    void volume_element_calc_geometry_properties(ElementType e_type,
                                                 const Index *element,
                                                 const Vector<Vec3> &nodes,
                                                 Scalar &volume,
                                                 Vec3 &centroid);

    void face_element_calc_centroid(ElementType e_type,
                                    const Index *element,
                                    const Vector<Vec3> &nodes,
                                    Vec3 &centroid);
    void face_element_calc_face_normal(ElementType e_type, const Index *element, const Vector<Vec3> &nodes, Vec3 &S_ij);
    void tetrahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                              Scalar &volume, Vec3 &centroid);
    void tetrahedron_calc_volume(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Scalar &volume);
    void tetrahedron_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid);

    void triangle_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &S_ij);
    void triangle_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &centroid);

    void pyramid_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, const Vec3 &n4,
                                          Scalar &volume, Vec3 &centroid);

    void hexahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                             const Vec3 &n4, const Vec3 &n5, const Vec3 &n6, const Vec3 &n7,
                                             Scalar &volume, Vec3 &centroid);
    void quadrilateral_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &S_ij);
    void quadrilateral_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid);

    /*A structure of arrays (SoA) containing the faces and their required properties*/
    class Faces
    {
        friend class FV_Grid;

        struct CellPair
        {
            CellPair(Index i, Index j) : i{i}, j{j} {} // For some reason I had to define the constructor to make emplace_back work
            CellPair() = default;
            Index i, j;
            bool operator<(CellPair other) const
            {
                if (i != other.i)
                    return i < other.i;
                assert(j != other.j); // Two faces can't have the same owner and neigbour cells
                return j < other.j;
            }
        };

        Vector<CellPair> cell_indices;
        Vector<Vec3> face_normals;
        Vector<Vec3> centroid_to_face_i;
        Vector<Vec3> centroid_to_face_j;

        void reserve(Index size);
        void resize_geometry_properties();
        void sort(Index begin, Index end, const Elements &face_elements_old, Elements &face_elements_to_sort);

    public:
        Index size() const { return cell_indices.size(); }

        Index get_cell_i(Index face_index) const { return cell_indices[face_index].i; }
        Index get_cell_j(Index face_index) const { return cell_indices[face_index].j; }
        const Vec3 &get_face_normal(Index face_index) const { return face_normals[face_index]; }
        const Vec3 &get_centroid_to_face_i(Index face_index) const { return centroid_to_face_i[face_index]; }
        const Vec3 &get_centroid_to_face_j(Index face_index) const { return centroid_to_face_j[face_index]; }
        string to_string(Index ij) const;
    };

    /*A structure of arrays (SoA) containing the cells and their required properties*/
    class Cells
    {
        friend class FV_Grid;
        Vector<Scalar> volumes;
        Vector<Vec3> centroids;

        void resize(Index size);
        void reserve(Index size);

    public:
        Index size() const { return volumes.size(); }
        Scalar get_cell_volume(Index cell_index) const { return volumes[cell_index]; }
        const Vec3 &get_centroid(Index cell_index) const { return centroids[cell_index]; }
        string to_string(Index i) const;
    };

    struct Patch
    {
        BoundaryType boundary_type;
        // Vector<Index> boundary_face_indices;
        Index N_FACES;
        Index FIRST_FACE;
    };

}