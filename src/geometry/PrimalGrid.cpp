#include "../../include/geometry/PrimalGrid.hpp"

namespace geometry
{

    void PrimalGrid::read_mesh(const Config &config)
    {
        string mesh_filename_path = config.get_mesh_filename_path();
        string extension = mesh_filename_path.substr(mesh_filename_path.find_last_of(".") + 1);
        if (extension == "nf")
            try
            {
                read_netflux_mesh(config);
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error("Error reading netflux mesh '" + mesh_filename_path + "', " + string(e.what()));
            }
        else if (extension == "su2")
            try
            {
                read_su2_mesh(config);
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error("Error reading su2 mesh '" + mesh_filename_path + "', " + string(e.what()));
            }
        else
            assert(false);

        assert(nodes.capacity() == nodes.size());
        volume_elements.shrink_to_fit();
        for (auto &element_patch : element_patches)
            element_patch.surface_elements.shrink_to_fit();
        assert(element_patches.capacity() == element_patches.size());

        cout << "Native mesh has been read\n";
    }

    void PrimalGrid::read_netflux_mesh(const Config &config)
    {
        assert(false); // Not updatet
        // string mesh_filename = config.get_mesh_filename_path();
        // std::ifstream ist{mesh_filename};
        // std::stringstream ss;

        // if (!ist)
        //     throw std::runtime_error("Couldn't open the mesh file: " + mesh_filename);

        // string line, tmp;
        // Index N_NODES, N_ELEMENTS;

        // // READ N_NODES
        // ist >> tmp >> N_NODES;
        // // FAIL_IF(tmp != "N_NODES");

        // // READ N_ELEMENTS
        // ist >> tmp;
        // ist >> N_ELEMENTS;
        // // FAIL_IF(tmp != "N_ELEMENTS");

        // nodes.resize(N_NODES);
        // while (getline(ist, line))
        // {
        //     if (line.size() != 0)
        //     {
        //         // FAIL_IF(line != "#nodes");
        //         break;
        //     }
        // }
        // // reading nodes
        // for (Index i{0}; i < N_NODES; i++)
        // {
        //     Vec3 node;
        //     ist >> node.x() >> node.y() >> node.z();
        //     nodes.at(i) = node;
        // }

        // tet_connect.resize(N_ELEMENTS);
        // while (getline(ist, line))
        // {
        //     if (line.size() != 0)
        //     {
        //         // FAIL_IF(line != "#tetrahedra");
        //         break;
        //     }
        // }
        // // reading elements (assuming only tetrahedra for now)
        // for (Index i{0}; i < N_ELEMENTS; i++)
        // {
        //     TetConnect t;
        //     ist >> t.a() >> t.b() >> t.c() >> t.d();
        //     tet_connect.at(i) = t;
        // }

        // while (getline(ist, line))
        // {
        //     if (line.size() != 0)
        //     {
        //         // FAIL_IF(line != "#patches");
        //         break;
        //     }
        // }

        // // Reading boundary patches
        // string patch_name;
        // Index N_surface_elements;
        // while (ist >> patch_name >> N_surface_elements)
        // {
        //     TriConnect tri;
        //     TriPatchConnect p;
        //     p.patch_name = patch_name;
        //     if (!config.input_file_contains_patch_name(p.patch_name))
        //         throw std::runtime_error("Patch with name '" + p.patch_name + "' is not named in input file\n");
        //     p.triangles.resize(N_surface_elements);
        //     for (Index i{0}; i < N_surface_elements; i++)
        //     {
        //         ist >> tri.a() >> tri.b() >> tri.c();
        //         p.triangles.at(i) = tri;
        //     }
        //     tri_patch_connect_list.emplace_back(p);
        // }
    }

    void PrimalGrid::read_su2_mesh(const Config &config)
    {
        const string mesh_filename = config.get_mesh_filename_path();
        std::ifstream ist{mesh_filename};
        std::stringstream ss;

        if (!ist)
            throw std::runtime_error("Couldn't open the mesh file: " + mesh_filename);

        string tmp_string;
        ShortIndex tmp_int, su2_e_id_int;
        Index N_NODES, N_ELEMENTS, N_PATCHES, element_num;

        auto check_string_correctness = [](string actual_string, string correct_string)
        {
            if (actual_string != correct_string)
                throw std::runtime_error("symbol '" + actual_string + "' parsed instead of correct symbol '" + correct_string + "'\n");
        };

        ist >> tmp_string >> tmp_int;
        check_string_correctness(tmp_string, "NDIME=");
        if (tmp_int != 3)
            throw std::runtime_error("Only three dimensions (NDIME=3) is accepted, NDIME = " + tmp_int);

        /*--------------------------------------------------------------------
        Reading element connectivity
        --------------------------------------------------------------------*/
        ist >> tmp_string >> N_ELEMENTS;
        check_string_correctness(tmp_string, "NELEM=");
        volume_elements.reserve(N_ELEMENTS, MAX_NODES_VOLUME_ELEMENT);
        array<Index, MAX_NODES_VOLUME_ELEMENT> volume_element_nodes;
        for (Index i{0}; i < N_ELEMENTS; i++)
        {
            ist >> su2_e_id_int;
            auto su2_e_id = static_cast<SU2_ElementID>(su2_e_id_int);
            if (su2_identifier_to_element_type.count(su2_e_id) == 0)
                throw std::runtime_error("Illegal su2 element identifier (" + std::to_string(su2_e_id_int) +
                                         ") detected in su2 mesh file\n");
            ElementType e_type = su2_identifier_to_element_type.at(su2_e_id);
            if (!is_volume_element.at(e_type))
                throw std::runtime_error("Face element of type " + element_type_to_string.at(e_type) +
                                         "detected in su2 mesh file in the volume mesh");
            for (ShortIndex k{0}; k < num_nodes_in_element.at(e_type); k++)
                ist >> volume_element_nodes[k];
            volume_elements.add_element(e_type, volume_element_nodes.data());
        }

        /*--------------------------------------------------------------------
        Reading nodes
        --------------------------------------------------------------------*/
        ist >> tmp_string >> N_NODES;
        check_string_correctness(tmp_string, "NPOIN=");
        nodes.resize(N_NODES);

        for (Index i{0}; i < N_NODES; i++)
        {
            Vec3 node;
            ist >> node.x() >> node.y() >> node.z();
            nodes.at(i) = node;
        }

        /*--------------------------------------------------------------------
        Reading boundary patches
        --------------------------------------------------------------------*/
        ist >> tmp_string >> N_PATCHES;
        check_string_correctness(tmp_string, "NMARK=");
        element_patches.resize(N_PATCHES);
        array<Index, MAX_NODES_FACE_ELEMENT> boundary_element_nodes;
        for (Index i{0}; i < N_PATCHES; i++)
        {
            Elements &boundary_elements = element_patches[i].boundary_elements;
            string &patch_name = element_patches[i].patch_name;

            ist >> tmp_string >> patch_name;
            check_string_correctness(tmp_string, "MARKER_TAG=");
            if (!config.input_file_contains_patch_name(patch_name))
                throw std::runtime_error("Patch with name '" + patch_name + "' is not named in input file\n");

            Index N_MARKER_ELEMENTS;
            ist >> tmp_string >> N_MARKER_ELEMENTS;
            check_string_correctness(tmp_string, "MARKER_ELEMS=");

            boundary_elements.reserve(N_MARKER_ELEMENTS, MAX_NODES_FACE_ELEMENT);
            for (Index j{0}; j < N_MARKER_ELEMENTS; j++)
            {
                ist >> su2_e_id_int;
                auto su2_e_id = static_cast<SU2_ElementID>(su2_e_id_int);
                if (su2_identifier_to_element_type.count(su2_e_id) == 0)
                    throw std::runtime_error("Illegal su2 element identifier (" + std::to_string(su2_e_id_int) +
                                             ") detected in su2 mesh file\n");
                ElementType e_type = su2_identifier_to_element_type.at(su2_e_id);
                if (is_volume_element.at(e_type))
                    throw std::runtime_error("Volume element of type " + element_type_to_string.at(e_type) +
                                             "detected in su2 mesh file among marker elements");
                for (ShortIndex k{0}; k < num_nodes_in_element.at(e_type); k++)
                    ist >> boundary_element_nodes[k];
                boundary_elements.add_element(e_type, boundary_element_nodes.data());
            }
        }
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

    inline FaceElement get_face_element_k_of_volume_element(ElementType volume_element_type,
                                                            const Index *volume_element,
                                                            ShortIndex face_k)
    {
        assert(volume_element_type == ElementType::Tetrahedron || volume_element_type == ElementType::Hexahedron);
        assert(is_volume_element.at(volume_element_type));
        assert(face_k <= num_nodes_in_element.at(volume_element_type));

        ElementType face_element_type;
        array<Index, MAX_NODES_FACE_ELEMENT> face_element;

        if (volume_element_type == ElementType::Tetrahedron)
            get_face_element_k_of_tetrahedron(volume_element_type, volume_element, face_k, face_element_type, face_element);
        else if (volume_element_type == ElementType::Hexahedron)
            get_face_element_k_of_hexahedron(volume_element_type, volume_element, face_k, face_element_type, face_element);
        return FaceElement{face_element_type, face_element.data()};
    }

    inline void get_face_element_k_of_tetrahedron(ElementType volume_element_type,
                                                  const Index *ve,
                                                  ShortIndex face_k,
                                                  ElementType &face_element_type,
                                                  array<Index, MAX_NODES_FACE_ELEMENT> &fe)
    {
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

    inline void get_face_element_k_of_hexahedron(ElementType volume_element_type,
                                                 const Index *ve,
                                                 ShortIndex face_k,
                                                 ElementType &face_element_type,
                                                 array<Index, MAX_NODES_FACE_ELEMENT> &fe)
    {
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

    inline void volume_element_calc_geometry_properties(ElementType e_type,
                                                        const Index *element,
                                                        const Vector<Vec3> &nodes,
                                                        Scalar &volume,
                                                        Vec3 &centroid)
    {
        assert(is_volume_element.at(e_type));
        if (e_type == ElementType::Tetrahedron)
            tetrahedron_calc_geometry_properties(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]],
                                                 volume, centroid);
        else if (e_type == ElementType::Hexahedron)
            hexahedron_calc_geometry_properties(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]],
                                                nodes[element[4]], nodes[element[5]], nodes[element[6]], nodes[element[7]],
                                                volume, centroid);
        assert(false); // invalid element type
    }

    inline void face_element_calc_centroid(ElementType e_type, const Index *element, const Vector<Vec3> &nodes, Vec3 &centroid)
    {
        assert(!is_volume_element.at(e_type));
        if (e_type == ElementType::Triangle)
            triangle_calc_centroid(nodes[element[0]], nodes[element[1]], nodes[element[2]], centroid);
        else if (e_type == ElementType::Quadrilateral)
            quadrilateral_calc_centroid(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]], centroid);
        assert(false); // Invalid element type
    }

    inline void face_element_calc_face_normal(ElementType e_type, const Index *element, const Vector<Vec3> &nodes, Vec3 &S_ij)
    {
        assert(!is_volume_element.at(e_type));
        if (e_type == ElementType::Triangle)
            triangle_calc_face_normal(nodes[element[0]], nodes[element[1]], nodes[element[2]], S_ij);
        else if (e_type == ElementType::Quadrilateral)
            quadrilateral_calc_face_normal(nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]], S_ij);
        assert(false); // Invalid element type
    }

    inline void tetrahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
                                                     Scalar &volume, Vec3 &centroid)
    {
        tetrahedron_calc_volume(n0, n1, n2, n3, volume);
        tetrahedron_calc_centroid(n0, n1, n2, n3, centroid);
    }

    inline void tetrahedron_calc_volume(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Scalar &volume)
    {
        Vec3 ab = n1 - n0;
        Vec3 ac = n2 - n0;
        Vec3 ad = n3 - n0;
        volume = -ab.cross(ac).dot(ad) / 6;
        assert(volume > 0);
    }
    inline void tetrahedron_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid)
    {
        centroid = 1.0 / 4 * (n0 + n1 + n2 + n3);
    }

    inline void triangle_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &S_ij)
    {
        /*Using that the cross product of two vectors is the area of a parallelogram (with correct normal),
        so the half is the area of a triangle*/
        Vec3 ab = n1 - n0;
        Vec3 ac = n2 - n0;
        S_ij = 0.5 * ab.cross(ac);
    }
    inline void triangle_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, Vec3 &centroid)
    {
        centroid = (n0 + n1 + n2) / 3;
    }

    inline void pyramid_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, const Vec3 &n4,
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

    inline void hexahedron_calc_geometry_properties(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3,
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

    inline void quadrilateral_calc_face_normal(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &S_ij)
    {
        Vec3 S_tri_a, S_tri_b;
        triangle_calc_face_normal(n0, n1, n2, S_tri_a);
        triangle_calc_face_normal(n0, n1, n3, S_tri_b);
        S_ij = 0.5 * (S_tri_a + S_tri_b);
    }

    inline void quadrilateral_calc_centroid(const Vec3 &n0, const Vec3 &n1, const Vec3 &n2, const Vec3 &n3, Vec3 &centroid)
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
}