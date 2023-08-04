#include "../../include/geometry/PrimalGrid.hpp"

namespace geometry
{

    PrimalGrid::PrimalGrid(const Config &config)
    {
        cerr << "PrimalGrid mesh read constructor called, rank " << NF_MPI::get_rank() << endl;
        try
        {
            read_mesh(config);
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error("Error creating primal grid:\n" + string(e.what()));
        }
    }

    PrimalGrid::PrimalGrid(Vector<Vec3> &&nodes, Elements &&vol_elements, Vector<ElementPatch> &&element_patches, Index eID_glob_first)
        : nodes{move(nodes)}, vol_elements{move(vol_elements)}, element_patches{move(element_patches)}, eID_glob_first{eID_glob_first}
    {
        cerr << "PrimalGrid existing element constructor called, rank " << NF_MPI::get_rank() << endl;
    }

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
        vol_elements.shrink_to_fit();
        for (auto &element_PatchExt : element_patches)
            element_PatchExt.boundary_elements.shrink_to_fit();
        assert(element_patches.capacity() == element_patches.size());

        if (config.check_grid_validity())
            partial_validity_check();

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
        // string PatchExt_name;
        // Index N_surface_elements;
        // while (ist >> PatchExt_name >> N_surface_elements)
        // {
        //     TriConnect tri;
        //     TriPatchExtConnect p;
        //     p.PatchExt_name = PatchExt_name;
        //     if (!config.input_file_contains_PatchExt_name(p.PatchExt_name))
        //         throw std::runtime_error("PatchExt with name '" + p.PatchExt_name + "' is not named in input file\n");
        //     p.triangles.resize(N_surface_elements);
        //     for (Index i{0}; i < N_surface_elements; i++)
        //     {
        //         ist >> tri.a() >> tri.b() >> tri.c();
        //         p.triangles.at(i) = tri;
        //     }
        //     tri_PatchExt_connect_list.emplace_back(p);
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
        ShortIndex tmp_int, vtk_e_id_int;
        Index N_NODES, N_ELEMENTS, N_patches;

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
        vol_elements.reserve(N_ELEMENTS, MAX_NODES_VOLUME_ELEMENT);
        Array<Index, MAX_NODES_VOLUME_ELEMENT> volume_element_nodes;
        for (Index i{0}; i < N_ELEMENTS; i++)
        {
            ist >> vtk_e_id_int;
            if (!legal_vtk_element_identifier(vtk_e_id_int))
                throw std::runtime_error("Illegal su2 volume element identifier (" + std::to_string(vtk_e_id_int) +
                                         ") detected in su2 mesh file\n");
            auto e_type = static_cast<ElementType>(vtk_e_id_int);

            if (!is_volume_element(e_type))
                throw std::runtime_error("Face element of type " + get_element_string(e_type) +
                                         "detected in su2 mesh file in the volume mesh");
            for (ShortIndex k{0}; k < get_num_nodes_in_element(e_type); k++)
                ist >> volume_element_nodes[k];
            ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to the next line
            vol_elements.add_element(e_type, volume_element_nodes.data());
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
        ist >> tmp_string >> N_patches;
        check_string_correctness(tmp_string, "NMARK=");
        element_patches.resize(N_patches);
        Array<Index, MAX_NODES_FACE_ELEMENT> boundary_element_nodes;
        for (Index i{0}; i < N_patches; i++)
        {
            Elements &boundary_elements = element_patches[i].boundary_elements;
            string &PatchExt_name = element_patches[i].PatchExt_name;

            ist >> tmp_string >> PatchExt_name;
            check_string_correctness(tmp_string, "MARKER_TAG=");
            if (!config.input_file_contains_PatchExt_name(PatchExt_name))
                throw std::runtime_error("PatchExt with name '" + PatchExt_name + "' is not named in input file\n");

            Index N_MARKER_ELEMENTS;
            ist >> tmp_string >> N_MARKER_ELEMENTS;
            check_string_correctness(tmp_string, "MARKER_ELEMS=");

            boundary_elements.reserve(N_MARKER_ELEMENTS, MAX_NODES_FACE_ELEMENT);
            for (Index j{0}; j < N_MARKER_ELEMENTS; j++)
            {
                ist >> vtk_e_id_int;
                if (!legal_vtk_element_identifier(vtk_e_id_int))
                    throw std::runtime_error("Illegal su2 element identifier (" + std::to_string(vtk_e_id_int) +
                                             ") detected in su2 mesh file\n");
                auto e_type = static_cast<ElementType>(vtk_e_id_int);
                if (is_volume_element(e_type))
                    throw std::runtime_error("Volume element of type " + get_element_string(e_type) +
                                             "detected in su2 mesh file among marker elements");
                for (ShortIndex k{0}; k < get_num_nodes_in_element(e_type); k++)
                    ist >> boundary_element_nodes[k];
                boundary_elements.add_element(e_type, boundary_element_nodes.data());
                ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to the next line
            }
        }
        /*Converting connectivities to vtk format*/
        vol_elements.salome_to_vtk_connectivity();
        for (auto &PatchExt : element_patches)
            PatchExt.boundary_elements.salome_to_vtk_connectivity();
    }
    void PrimalGrid::print_grid() const

    {
        cout << "NODES:\n";
        for (Index i{0}; i < nodes.size(); i++)
        {
            cout << i << ": " << array_to_string(nodes.at(i).data(), N_DIM) << endl;
        }
        cout << "\n\nVOLUME ELEMENTS:\n";
        for (Index i{0}; i < vol_elements.size(); i++)
        {
            cout << vol_elements.to_string(i);
        }
        cout << "\n\nPatchExt ELEMENTS:\n";
        for (const auto &ep : element_patches)
        {
            cout << "BC type: " << ep.PatchExt_name << endl;
            for (Index i{0}; i < ep.boundary_elements.size(); i++)
            {
                cout << i << ": " << ep.boundary_elements.to_string(i) << endl;
            }
            cout << "\n\n";
        }
    }
    /*Checking if volume is positive and if the centroid is within bounding box of element nodes*/
    void PrimalGrid::partial_validity_check()
    {
        cout << "Running some checks on primal mesh..\n";
        for (Index i{0}; i < vol_elements.size(); i++)
        {
            Vec3 centroid;
            Scalar vol;
            // cerr << "elem " << get_element_string(vol_elements.get_element_type(i))
            //      << " " << array_to_string(vol_elements.get_element_nodes(i), vol_elements.get_n_element_nodes(i)) << endl;

            volume_element_calc_geometry_properties(vol_elements.get_element_type(i),
                                                    vol_elements.get_element_nodes(i),
                                                    nodes, vol, centroid);

            Vec3 max_coord, min_coord;
            max_coord.setConstant(-std::numeric_limits<double>::infinity());
            min_coord.setConstant(std::numeric_limits<double>::infinity());

            for (ShortIndex k{0}; k < vol_elements.get_n_element_nodes(i); k++)
            {
                const Vec3 &node = nodes[vol_elements.get_element_nodes(i)[k]];
                max_coord = max_coord.cwiseMax(centroid.cwiseMax(node));
                min_coord = min_coord.cwiseMin(centroid.cwiseMin(node));
            }

            if (vol < 0.0)
                std::runtime_error("Element with negative volume detected in primal mesh. (Element " + std::to_string(i) + ")\n");

            for (Index i_dim{0}; i_dim < N_DIM; i_dim++)
                if (is_approx_equal(max_coord[i_dim], centroid[i_dim]) || is_approx_equal(max_coord[i_dim], centroid[i_dim]))
                    std::runtime_error("Element with centroid outside bounding box detected in primal mesh. (Element " + std::to_string(i) + ")\n");
        }
    }
    Index PrimalGrid::find_num_ghost_external() const
    {
        Index num_ghost{0};
        for (const auto &ep : element_patches)
            num_ghost += ep.boundary_elements.size();
        return num_ghost;
    }
}