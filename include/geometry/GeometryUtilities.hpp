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
        Pyramid = 14,
        Wedge = 13
    };

    constexpr ShortIndex N_NODES_TET{4};
    constexpr ShortIndex N_FACES_TET{4};

    constexpr ShortIndex N_NODES_TRI{3};

    constexpr ShortIndex N_NODES_HEX{8};
    constexpr ShortIndex N_FACES_HEX{6};

    constexpr ShortIndex N_NODES_QUAD{4};

    constexpr ShortIndex N_NODES_PYRAMID{5};
    constexpr ShortIndex N_FACES_PYRAMID{5};

    constexpr ShortIndex N_NODES_WEDGE{6};
    constexpr ShortIndex N_FACES_WEDGE{5};

    /*These constants store the maximum number of nodes among all element types
    for volume and face elements*/
    constexpr ShortIndex MAX_NODES_VOLUME_ELEMENT{N_NODES_HEX};
    constexpr ShortIndex MAX_NODES_FACE_ELEMENT{N_NODES_QUAD};
    namespace element_details
    {
        const map<ElementType, ShortIndex> element_num_nodes = {{ElementType::Triangle, N_NODES_TRI},
                                                                {ElementType::Quadrilateral, N_NODES_QUAD},
                                                                {ElementType::Tetrahedron, N_NODES_TET},
                                                                {ElementType::Hexahedron, N_NODES_HEX},
                                                                {ElementType::Pyramid, N_NODES_PYRAMID},
                                                                {ElementType::Wedge, N_NODES_WEDGE}};
        const map<ElementType, ShortIndex> volume_element_num_faces = {{ElementType::Tetrahedron, N_FACES_TET},
                                                                       {ElementType::Hexahedron, N_FACES_HEX},
                                                                       {ElementType::Pyramid, N_FACES_PYRAMID},
                                                                       {ElementType::Wedge, N_FACES_WEDGE}};
        const map<ElementType, string> element_to_string = {{ElementType::Triangle, "Triangle"},
                                                            {ElementType::Quadrilateral, "Quadrilateral"},
                                                            {ElementType::Tetrahedron, "Tetrahedron"},
                                                            {ElementType::Hexahedron, "Hexahedron"},
                                                            {ElementType::Pyramid, "Pyramid"},
                                                            {ElementType::Wedge, "Wedge"}};
        const map<ElementType, bool> element_is_volume = {{ElementType::Triangle, false},
                                                          {ElementType::Tetrahedron, true},
                                                          {ElementType::Quadrilateral, false},
                                                          {ElementType::Hexahedron, true},
                                                          {ElementType::Pyramid, true},
                                                          {ElementType::Wedge, true}};
    }
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
        case ElementType::Wedge:
            return true;
        default:
            return false;
        }
    }
    inline ShortIndex get_num_nodes_in_element(ElementType e_type)
    {
        assert(element_details::element_num_nodes.count(e_type) == 1);
        return element_details::element_num_nodes.at(e_type);
    }

    inline string get_element_string(ElementType e_type)
    {
        assert(element_details::element_to_string.count(e_type) == 1);
        return element_details::element_to_string.at(e_type);
    }

    inline bool is_volume_element(ElementType e_type)
    {
        assert(element_details::element_is_volume.count(e_type) == 1);
        return element_details::element_is_volume.at(e_type);
    }
    inline ShortIndex get_num_faces_volume_element(ElementType e_type)
    {
        assert(is_volume_element(e_type));
        assert(element_details::volume_element_num_faces.count(e_type) == 1);
        return element_details::volume_element_num_faces.at(e_type);
    }

    /*Connectivities of elements laid out on the compressed row storage format (CRS) as METIS uses.*/
    class Elements
    {
        friend class Faces;

    protected:
        Vector<Index> e_ptr = {0};
        Vector<Index> e_ind;
        Vector<ElementType> element_types;

    public:
        Elements() = default;
        Elements(const Elements &other);
        Elements &operator=(Elements other);
        Index size() const { return e_ptr.size() - 1; }

        const Vector<Index> &get_e_ptr() const { return e_ptr; }
        const Vector<Index> &get_e_ind() const { return e_ind; }
        Vector<Index> &get_e_ptr() { return e_ptr; }
        Vector<Index> &get_e_ind() { return e_ind; }

        void add_element(ElementType e_type, const Index *element)
        {
            ShortIndex n_nodes = get_num_nodes_in_element(e_type);
            e_ptr.emplace_back(e_ptr.back() + n_nodes);
            for (ShortIndex k{0}; k < n_nodes; k++)
                e_ind.emplace_back(element[k]);
            element_types.emplace_back(e_type);
        }
        /*Special function that uses a map to convert from global to local node indices
        when using multiple processors*/
        void add_element_local(ElementType e_type, const Index *element, const map<Index, Index> &glob_to_loc)
        {
            ShortIndex n_nodes = get_num_nodes_in_element(e_type);
            e_ptr.emplace_back(e_ptr.back() + n_nodes);
            for (ShortIndex k{0}; k < n_nodes; k++)
            {
                assert(glob_to_loc.count(element[k] == 1));
                e_ind.emplace_back(glob_to_loc.at(element[k]));
            }
            element_types.emplace_back(e_type);
        }

        Index sum_nodes_over_all_elements() const { return e_ptr.back(); }

        const Index *get_element_nodes(Index i) const
        {
            assert((e_ptr[i + 1] - e_ptr[i]) == get_num_nodes_in_element(element_types[i]));
            return &e_ind[e_ptr[i]];
        }

    private:
        Index *get_element_nodes_modifiable(Index i)
        {
            assert((e_ptr[i + 1] - e_ptr[i]) == get_num_nodes_in_element(element_types[i]));
            return &e_ind[e_ptr[i]];
        }

    public:
        const ShortIndex get_n_element_nodes(Index i) const
        {
            Index num_nodes{e_ptr[i + 1] - e_ptr[i]};
            assert(num_nodes == get_num_nodes_in_element(element_types[i]));
            return num_nodes;
        }

        ElementType get_element_type(Index i) const
        {
            return element_types[i];
        }

        void reserve(Index n_elements, ShortIndex max_nodes_per_element);
        void shrink_to_fit();

        string to_string(Index i) const;

        /*--------------------------------------------------------------------
        It seems like the mesh generated by salome and vtk uses different element
        connectivity definitions.
        Salome: https://docs.salome-platform.org/latest/gui/SMESH/connectivity.html
        vtk: http://www.princeton.edu/~efeibush/viscourse/vtk.pdf
        Since the vtk format is currently used, I add a function that reorders the
        elements from salome to vtk.
        --------------------------------------------------------------------*/
        void salome_to_vtk_connectivity();
    };

    class FaceElement
    {
        std::array<Index, MAX_NODES_FACE_ELEMENT> nodes_sorted;

    public:
        std::array<Index, MAX_NODES_FACE_ELEMENT> nodes;
        const ShortIndex n_nodes;
        const ElementType e_type;

        FaceElement(ElementType e_type, const Index *element) : n_nodes{get_num_nodes_in_element(e_type)}, e_type{e_type}

        {
            assert(!is_volume_element(e_type));
            std::copy(element, element + n_nodes, nodes.begin());
            std::copy(element, element + n_nodes, nodes_sorted.begin());
            std::sort(nodes_sorted.begin(), nodes_sorted.begin() + n_nodes);
        }

        bool operator<(const FaceElement &other) const
        {
            if (n_nodes != other.n_nodes)
                return n_nodes < other.n_nodes;

            for (ShortIndex i{0}; i < n_nodes; i++)
                if (nodes_sorted[i] != other.nodes_sorted[i])
                    return nodes_sorted[i] < other.nodes_sorted[i];
            return false;
        }
    };

    FaceElement get_face_element_k_of_volume_element(ElementType volume_element_type,
                                                     const Index *volume_element,
                                                     ShortIndex face_k);
    void get_face_element_k_of_tetrahedron(const Index *ve,
                                           ShortIndex face_k,
                                           ElementType &face_element_type,
                                           array<Index, MAX_NODES_FACE_ELEMENT> &fe);
    void get_face_element_k_of_hexahedron(const Index *ve,
                                          ShortIndex face_k,
                                          ElementType &face_element_type,
                                          array<Index, MAX_NODES_FACE_ELEMENT> &fe);
    void get_face_element_k_of_pyramid(const Index *ve,
                                       ShortIndex face_k,
                                       ElementType &face_element_type,
                                       array<Index, MAX_NODES_FACE_ELEMENT> &fe);
    void get_face_element_k_of_wedge(const Index *ve,
                                     ShortIndex face_k,
                                     ElementType &face_element_type,
                                     array<Index, MAX_NODES_FACE_ELEMENT> &fe);

    // struct SortedFaceElement : public FaceElement
    // {
    //     SortedFaceElement(const FaceElement &rhs) : FaceElement(rhs)
    //     {
    //         std::sort(nodes.begin(), nodes.begin() + n_nodes);
    //     }

    //     bool operator<(const SortedFaceElement &other) const
    //     {
    //         if (n_nodes == other.n_nodes)
    //         {
    //             for (ShortIndex i{0}; i < n_nodes; i++)
    //                 if (nodes[i] != other.nodes[i])
    //                     return nodes[i] < other.nodes[i];
    //             return false;
    //         }
    //         else
    //             return n_nodes < other.n_nodes;
    //     }
    // };

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

    void hexahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                             const Vec3 &n4, const Vec3 &n5, const Vec3 &n6, const Vec3 &n7,
                                             Scalar &volume, Vec3 &centroid);
    void pyramid_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, const Vec3 &n4,
                                          Scalar &volume, Vec3 &centroid);
    void wedge_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, const Vec3 &n4, const Vec3 &n5,
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