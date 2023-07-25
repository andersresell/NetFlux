#include "../../include/geometry/GeometryUtilities.hpp"

namespace geometry
{

    Elements::Elements(const Elements &other)

    {
        n_ptr.resize(other.n_ptr.size());
        n_ind.resize(other.n_ind.size());
        element_types.resize(other.element_types.size());
        std::copy(other.n_ptr.begin(), other.n_ptr.end(), n_ptr.begin());
        std::copy(other.n_ind.begin(), other.n_ind.end(), n_ind.begin());
        std::copy(other.element_types.begin(), other.element_types.end(), element_types.begin());
    }

    Elements &Elements::operator=(Elements other)
    {
        std::swap(this->n_ptr, other.n_ptr);
        std::swap(this->n_ind, other.n_ind);
        std::swap(this->element_types, other.element_types);
        return *this;
    }

    void Elements::reserve(Index n_elements, ShortIndex max_nodes_per_element)
    {
        n_ptr.reserve(n_elements + 1);
        n_ind.reserve(n_elements * max_nodes_per_element); // This will give some redundency if multiple element types are present in the mesh
        element_types.reserve(n_elements);
    }
    void Elements::shrink_to_fit()
    {
        n_ptr.shrink_to_fit();
        n_ind.shrink_to_fit();
        element_types.shrink_to_fit();
    }

    string Elements::to_string(Index i) const
    {
        ShortIndex n_nodes = get_n_element_nodes(i);
        const Index *nodes = get_element_nodes(i);
        std::stringstream ss;
        ss << array_to_string(nodes, n_nodes) << ", " << get_element_string(get_element_type(i)) << endl;
        return ss.str();
    }

    FaceElement get_face_element_k_of_volume_element(ElementType volume_element_type,
                                                     const Index *volume_element,
                                                     ShortIndex face_k)
    {
        assert(volume_element_type == ElementType::Tetrahedron || volume_element_type == ElementType::Hexahedron);
        assert(is_volume_element(volume_element_type));
        assert(face_k <= get_num_faces_volume_element(volume_element_type));
        ElementType face_element_type;
        array<Index, MAX_NODES_FACE_ELEMENT> face_element;

        if (volume_element_type == ElementType::Tetrahedron)
            get_face_element_k_of_tetrahedron(volume_element, face_k, face_element_type, face_element);
        else if (volume_element_type == ElementType::Hexahedron)
            get_face_element_k_of_hexahedron(volume_element, face_k, face_element_type, face_element);
        else if (volume_element_type == ElementType::Pyramid)
            get_face_element_k_of_pyramid(volume_element, face_k, face_element_type, face_element);
        else if (volume_element_type == ElementType::Wedge)
            get_face_element_k_of_wedge(volume_element, face_k, face_element_type, face_element);

        return FaceElement{face_element_type, face_element.data()};
    }

    void get_face_element_k_of_tetrahedron(const Index *ve,
                                           ShortIndex face_k,
                                           ElementType &face_element_type,
                                           array<Index, MAX_NODES_FACE_ELEMENT> &fe)
    {
        assert(face_k < N_FACES_TET);
        face_element_type = ElementType::Triangle;
        if (face_k == 0)
        {
            fe[0] = ve[0];
            fe[1] = ve[1];
            fe[2] = ve[2];
        }
        else if (face_k == 1)
        {
            fe[0] = ve[0];
            fe[1] = ve[1];
            fe[2] = ve[3];
        }
        else if (face_k == 2)
        {
            fe[0] = ve[0];
            fe[1] = ve[2];
            fe[2] = ve[3];
        }
        else if (face_k == 3)
        {
            fe[0] = ve[1];
            fe[1] = ve[2];
            fe[2] = ve[3];
        }
        /*When time: Check for possible branching in assembly output and consider branchless implementation (example below).
        (I suspect the compiler optimizes this away here)*/
        // fe[0] = ve[0]*(face_k==0) + ve[0]*(face_k==1) + ve[0]*(face_k==2)+ ve[3]*(face_k==3);
    }

    void get_face_element_k_of_hexahedron(const Index *ve,
                                          ShortIndex face_k,
                                          ElementType &face_element_type,
                                          array<Index, MAX_NODES_FACE_ELEMENT> &fe)
    {
        assert(face_k < N_FACES_HEX);
        face_element_type = ElementType::Quadrilateral;
        if (face_k == 0)
        {
            fe[0] = ve[0];
            fe[1] = ve[1];
            fe[2] = ve[2];
            fe[3] = ve[3];
        }
        else if (face_k == 1)
        {
            fe[0] = ve[4];
            fe[1] = ve[5];
            fe[2] = ve[6];
            fe[3] = ve[7];
        }
        else if (face_k == 2)
        {
            fe[0] = ve[0];
            fe[1] = ve[1];
            fe[2] = ve[5];
            fe[3] = ve[4];
        }
        else if (face_k == 3)
        {
            fe[0] = ve[3];
            fe[1] = ve[2];
            fe[2] = ve[6];
            fe[3] = ve[7];
        }
        else if (face_k == 4)
        {
            fe[0] = ve[3];
            fe[1] = ve[0];
            fe[2] = ve[4];
            fe[3] = ve[7];
        }
        else if (face_k == 5)
        {
            fe[0] = ve[1];
            fe[1] = ve[2];
            fe[2] = ve[6];
            fe[3] = ve[5];
        }
    }
    void get_face_element_k_of_pyramid(const Index *ve,
                                       ShortIndex face_k,
                                       ElementType &face_element_type,
                                       array<Index, MAX_NODES_FACE_ELEMENT> &fe)
    {
        assert(face_k < N_FACES_PYRAMID);
        face_element_type = (face_k == 0) ? ElementType::Quadrilateral : ElementType::Triangle;
        if (face_k == 0)
        {
            fe[0] = ve[0];
            fe[1] = ve[1];
            fe[2] = ve[2];
            fe[3] = ve[3];
        }
        else if (face_k == 1)
        {
            fe[0] = ve[0];
            fe[1] = ve[1];
            fe[2] = ve[4];
        }
        else if (face_k == 2)
        {
            fe[0] = ve[1];
            fe[1] = ve[2];
            fe[2] = ve[4];
        }
        else if (face_k == 3)
        {
            fe[0] = ve[2];
            fe[1] = ve[3];
            fe[2] = ve[4];
        }
        else if (face_k == 4)
        {
            fe[0] = ve[3];
            fe[1] = ve[0];
            fe[2] = ve[4];
        }
    }
    void get_face_element_k_of_wedge(const Index *ve,
                                     ShortIndex face_k,
                                     ElementType &face_element_type,
                                     array<Index, MAX_NODES_FACE_ELEMENT> &fe)
    {
        assert(face_k < N_FACES_WEDGE);
        if (face_k == 0)
        {
            face_element_type = ElementType::Triangle;
            fe[0] = ve[0];
            fe[1] = ve[1];
            fe[2] = ve[2];
        }
        else if (face_k == 1)
        {
            face_element_type = ElementType::Triangle;
            fe[0] = ve[4];
            fe[1] = ve[3];
            fe[2] = ve[5];
        }
        else if (face_k == 2)
        {
            face_element_type = ElementType::Quadrilateral;
            fe[0] = ve[1];
            fe[1] = ve[0];
            fe[2] = ve[3];
            fe[3] = ve[4];
        }
        else if (face_k == 3)
        {
            face_element_type = ElementType::Quadrilateral;
            fe[0] = ve[0];
            fe[1] = ve[2];
            fe[2] = ve[5];
            fe[3] = ve[3];
        }
        else if (face_k == 4)
        {
            face_element_type = ElementType::Quadrilateral;
            fe[0] = ve[4];
            fe[1] = ve[5];
            fe[2] = ve[2];
            fe[3] = ve[1];
        }
    }

    void volume_element_calc_geometry_properties(ElementType e_type,
                                                 const Index *element,
                                                 const Vector<Vec3> &nodes,
                                                 Scalar &volume,
                                                 Vec3 &centroid)
    {
        assert(is_volume_element(e_type));
        if (e_type == ElementType::Tetrahedron)
            tetrahedron_calc_geometry_properties(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]],
                                                 volume, centroid);
        else if (e_type == ElementType::Hexahedron)
            hexahedron_calc_geometry_properties(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]],
                                                nodes[element[4]], nodes[element[5]], nodes[element[6]], nodes[element[7]],
                                                volume, centroid);
        else if (e_type == ElementType::Pyramid)
            pyramid_calc_geometry_properties(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]],
                                             nodes[element[4]], volume, centroid);
        else if (e_type == ElementType::Wedge)
            wedge_calc_geometry_properties(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]],
                                           nodes[element[4]], nodes[element[5]], volume, centroid);
    }

    void face_element_calc_centroid(ElementType e_type, const Index *element, const Vector<Vec3> &nodes, Vec3 &centroid)
    {
        assert(e_type == ElementType::Triangle || e_type == ElementType::Quadrilateral);
        if (e_type == ElementType::Triangle)
            triangle_calc_centroid(nodes[element[0]], nodes[element[1]], nodes[element[2]], centroid);
        else if (e_type == ElementType::Quadrilateral)
            quadrilateral_calc_centroid(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]], centroid);
    }

    void face_element_calc_face_normal(ElementType e_type, const Index *element, const Vector<Vec3> &nodes, Vec3 &S_ij)
    {
        assert(e_type == ElementType::Triangle || e_type == ElementType::Quadrilateral);
        if (e_type == ElementType::Triangle)
            triangle_calc_face_normal(nodes[element[0]], nodes[element[1]], nodes[element[2]], S_ij);
        else if (e_type == ElementType::Quadrilateral)
            quadrilateral_calc_face_normal(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]], S_ij);
    }

    void tetrahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                              Scalar &volume, Vec3 &centroid)
    {
        tetrahedron_calc_volume(n0, n1, n2, n3, volume);
        tetrahedron_calc_centroid(n0, n1, n2, n3, centroid);
    }

    void tetrahedron_calc_volume(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Scalar &volume)
    {
        Vec3 ab = n1 - n0;
        Vec3 ac = n2 - n0;
        Vec3 ad = n3 - n0;
        volume = -ab.cross(ac).dot(ad) / 6;
        assert(volume > 0);
    }
    void tetrahedron_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid)
    {
        centroid = 1.0 / 4 * (n0 + n1 + n2 + n3);
    }

    void triangle_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &S_ij)
    {
        /*Using that the cross product of two vectors is the area of a parallelogram (with correct normal),
        so the half is the area of a triangle*/
        Vec3 ab = n1 - n0;
        Vec3 ac = n2 - n0;
        S_ij = 0.5 * ab.cross(ac);
    }
    void triangle_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &centroid)
    {
        centroid = (n0 + n1 + n2) / 3;
    }

    void pyramid_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, const Vec3 &n4,
                                          Scalar &volume, Vec3 &centroid)
    {
        Vec3 centroid_base;
        quadrilateral_calc_centroid(n0, n1, n2, n3, centroid_base);
        Vec3 distance_top_to_centroid_base = centroid_base - n4;
        Vec3 S_base;
        quadrilateral_calc_face_normal(n0, n1, n2, n3, S_base);

        centroid = 0.75 * centroid_base + 0.25 * n4;
        volume = distance_top_to_centroid_base.dot(S_base) / 3;
        assert(volume > 0);
    }
    void wedge_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, const Vec3 &n4, const Vec3 &n5,
                                        Scalar &volume, Vec3 &centroid)
    {
        volume = 0;
        centroid = {0, 0, 0};
        Vec3 centroid_geometric = (n0 + n1 + n2 + n3 + n4 + n5) / 6;
        Scalar vol_subgeom;
        Vec3 centroid_subgeom;

        tetrahedron_calc_geometry_properties(n0, n2, n1, centroid_geometric, vol_subgeom, centroid_subgeom);
        volume += vol_subgeom;
        centroid += vol_subgeom * centroid_subgeom;

        tetrahedron_calc_geometry_properties(n3, n4, n5, centroid_geometric, vol_subgeom, centroid_subgeom);
        volume += vol_subgeom;
        centroid += vol_subgeom * centroid_subgeom;

        pyramid_calc_geometry_properties(n0, n3, n5, n2, centroid_geometric, vol_subgeom, centroid_subgeom);
        volume += vol_subgeom;
        centroid += vol_subgeom * centroid_subgeom;

        pyramid_calc_geometry_properties(n0, n1, n4, n3, centroid_geometric, vol_subgeom, centroid_subgeom);
        volume += vol_subgeom;
        centroid += vol_subgeom * centroid_subgeom;

        pyramid_calc_geometry_properties(n0, n1, n5, n4, centroid_geometric, vol_subgeom, centroid_subgeom);
        volume += vol_subgeom;
        centroid += vol_subgeom * centroid_subgeom;

        centroid /= volume;
    }

    void hexahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                             const Vec3 &n4, const Vec3 &n5, const Vec3 &n6, const Vec3 &n7,
                                             Scalar &volume, Vec3 &centroid)
    {
        /*Composing the hexahedron into 6 pyramids and summing up each individual volume (Following Moukalled book)*/
        volume = 0;
        centroid = {0, 0, 0};
        Vec3 centroid_geometric = (n0 + n1 + n2 + n3 + n4 + n5 + n6 + n7) / 8; // Geometric centre of quad
        Scalar volume_pyramid;
        Vec3 centroid_pyramid;

        pyramid_calc_geometry_properties(n0, n1, n2, n3, centroid_geometric, volume_pyramid, centroid_pyramid);
        volume += volume_pyramid;
        centroid += volume_pyramid * centroid_pyramid;

        pyramid_calc_geometry_properties(n4, n5, n6, n7, centroid_geometric, volume_pyramid, centroid_pyramid);
        volume += volume_pyramid;
        centroid += volume_pyramid * centroid_pyramid;

        pyramid_calc_geometry_properties(n3, n0, n4, n7, centroid_geometric, volume_pyramid, centroid_pyramid);
        volume += volume_pyramid;
        centroid += volume_pyramid * centroid_pyramid;

        pyramid_calc_geometry_properties(n1, n2, n6, n5, centroid_geometric, volume_pyramid, centroid_pyramid);
        volume += volume_pyramid;
        centroid += volume_pyramid * centroid_pyramid;

        pyramid_calc_geometry_properties(n0, n1, n5, n4, centroid_geometric, volume_pyramid, centroid_pyramid);
        volume += volume_pyramid;
        centroid += volume_pyramid * centroid_pyramid;

        pyramid_calc_geometry_properties(n2, n3, n7, n6, centroid_geometric, volume_pyramid, centroid_pyramid);
        volume += volume_pyramid;
        centroid += volume_pyramid * centroid_pyramid;

        centroid /= volume;
    }

    void quadrilateral_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &S_ij)
    {
        Vec3 S_tri_a, S_tri_b;
        triangle_calc_face_normal(n0, n1, n2, S_tri_a);
        triangle_calc_face_normal(n0, n1, n3, S_tri_b);
        S_ij = 0.5 * (S_tri_a + S_tri_b);
    }

    void quadrilateral_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid)
    {
        Vec3 S_ij_tri_a, S_ij_tri_b;
        triangle_calc_face_normal(n0, n1, n2, S_ij_tri_a);
        triangle_calc_face_normal(n0, n2, n3, S_ij_tri_b);

        Vec3 centroid_geometric_tri_a, centroid_geometric_tri_b;
        triangle_calc_centroid(n0, n1, n2, centroid_geometric_tri_a);
        triangle_calc_centroid(n0, n2, n3, centroid_geometric_tri_b);

        Scalar S_tri_a = S_ij_tri_a.norm();
        Scalar S_tri_b = S_ij_tri_b.norm();
        Scalar S_quad = S_tri_a + S_tri_b;
        centroid = (S_tri_a * centroid_geometric_tri_a + S_tri_b * centroid_geometric_tri_b) / S_quad;
    }

    void Faces::reserve(Index size)
    {
        cell_indices.reserve(size);
        face_normals.reserve(size);
        centroid_to_face_i.reserve(size);
        centroid_to_face_j.reserve(size);
    }

    void Faces::resize_geometry_properties()
    {
        face_normals.resize(cell_indices.size());
        centroid_to_face_i.resize(cell_indices.size());
        centroid_to_face_j.resize(cell_indices.size());
    }

    /*Sorting all the Vectors from indices begin to end based on the cell_indices*/
    void Faces::sort(Index begin, Index end, const Elements &face_elements_old, Elements &face_elements_to_sort)
    {
        assert(begin < end && end <= cell_indices.size());
        assert(face_elements_to_sort.size() == begin);

        // Vector<Index> indices(end - begin);
        // for (Index i{begin}; i < end; i++)
        //     indices[i - begin] = i;

        Vector<Index> indices(end - begin);
        for (Index i{begin}; i < end; i++)
            indices[i - begin] = i;

        std::sort(indices.begin(), indices.end(), [this](Index a, Index b)
                  { return cell_indices[a] < cell_indices[b]; });

        /*Copying all content of faces */
        Vector<CellPair> cell_indices_copy(end - begin);
        Vector<Vec3> face_normals_copy(end - begin);
        Vector<Vec3> centroid_to_face_i_copy(end - begin);
        Vector<Vec3> centroid_to_face_j_copy(end - begin);
        std::copy(cell_indices.begin() + begin, cell_indices.begin() + end, cell_indices_copy.begin());
        std::copy(face_normals.begin() + begin, face_normals.begin() + end, face_normals_copy.begin());
        std::copy(centroid_to_face_i.begin() + begin, centroid_to_face_i.begin() + end, centroid_to_face_i_copy.begin());
        std::copy(centroid_to_face_j.begin() + begin, centroid_to_face_j.begin() + end, centroid_to_face_j_copy.begin());

        for (Index i{begin}; i < end; i++)
        {
            cell_indices[i] = cell_indices_copy[indices[i - begin] - begin];
            face_normals[i] = face_normals_copy[indices[i - begin] - begin];
            centroid_to_face_i[i] = centroid_to_face_i_copy[indices[i - begin] - begin];
            centroid_to_face_j[i] = centroid_to_face_j_copy[indices[i - begin] - begin];
            face_elements_to_sort.add_element(face_elements_old.get_element_type(indices[i - begin]),
                                              face_elements_old.get_element_nodes(indices[i - begin]));
        }
        assert(face_elements_to_sort.size() == end);
    }
    string Faces::to_string(Index ij) const
    {
        stringstream ss;
        ss << ij << ": "
           << "(" << cell_indices[ij].i << ", " << cell_indices[ij].j << "), S_ij = " << array_to_string(face_normals[ij].data(), N_DIM);
        return ss.str();
    }

    void Cells::resize(Index size)
    {
        volumes.resize(size);
        centroids.resize(size);
    }
    void Cells::reserve(Index size)
    {
        volumes.reserve(size);
        centroids.reserve(size);
    }
    string Cells::to_string(Index i) const
    {
        stringstream ss;
        ss << i << ": vol = " << volumes[i] << ", centroid = " << array_to_string(centroids[i].data(), N_DIM);
        return ss.str();
    }
}