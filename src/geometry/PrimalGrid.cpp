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

        cout << "Native mesh has been read\n";
    }

    void PrimalGrid::read_netflux_mesh(const Config &config)
    {
        assert(false); // Fix error handling
        string mesh_filename = config.get_mesh_filename_path();
        std::ifstream ist{mesh_filename};
        std::stringstream ss;

        if (!ist)
            throw std::runtime_error("Couldn't open the mesh file: " + mesh_filename);

        string line, tmp;
        Index N_NODES, N_ELEMENTS;

        // READ N_NODES
        ist >> tmp >> N_NODES;
        // FAIL_IF(tmp != "N_NODES");

        // READ N_ELEMENTS
        ist >> tmp;
        ist >> N_ELEMENTS;
        // FAIL_IF(tmp != "N_ELEMENTS");

        nodes.resize(N_NODES);
        while (getline(ist, line))
        {
            if (line.size() != 0)
            {
                // FAIL_IF(line != "#nodes");
                break;
            }
        }
        // reading nodes
        for (Index i{0}; i < N_NODES; i++)
        {
            Vec3 node;
            ist >> node.x() >> node.y() >> node.z();
            nodes.at(i) = node;
        }

        tet_connect.resize(N_ELEMENTS);
        while (getline(ist, line))
        {
            if (line.size() != 0)
            {
                // FAIL_IF(line != "#tetrahedra");
                break;
            }
        }
        // reading elements (assuming only tetrahedra for now)
        for (Index i{0}; i < N_ELEMENTS; i++)
        {
            TetConnect t;
            ist >> t.a() >> t.b() >> t.c() >> t.d();
            tet_connect.at(i) = t;
        }

        while (getline(ist, line))
        {
            if (line.size() != 0)
            {
                // FAIL_IF(line != "#patches");
                break;
            }
        }

        // Reading boundary patches
        string patch_name;
        Index N_surface_elements;
        while (ist >> patch_name >> N_surface_elements)
        {
            TriConnect tri;
            TriPatchConnect p;
            p.patch_name = patch_name;
            if (!config.input_file_contains_patch_name(p.patch_name))
                throw std::runtime_error("Patch with name '" + p.patch_name + "' is not named in input file\n");
            p.triangles.resize(N_surface_elements);
            for (Index i{0}; i < N_surface_elements; i++)
            {
                ist >> tri.a() >> tri.b() >> tri.c();
                p.triangles.at(i) = tri;
            }
            tri_patch_connect_list.emplace_back(p);
        }
    }

    void PrimalGrid::read_su2_mesh(const Config &config)
    {
        const string mesh_filename = config.get_mesh_filename_path();
        std::ifstream ist{mesh_filename};
        std::stringstream ss;

        if (!ist)
            throw std::runtime_error("Couldn't open the mesh file: " + mesh_filename);

        const ShortIndex SU2_TET_TYPE = 10, SU2_TRI_TYPE = 5;
        string tmp_string;
        ShortIndex tmp_int, element_type;
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
        // Reading element connectivity. Only permitting tetrahedral type
        ist >> tmp_string >> N_ELEMENTS;
        check_string_correctness(tmp_string, "NELEM=");
        tet_connect.resize(N_ELEMENTS);

        for (Index i{0}; i < N_ELEMENTS; i++)
        {
            TetConnect t;
            ist >> element_type >> t.a() >> t.b() >> t.c() >> t.d() >> element_num;
            if (element_type != SU2_TET_TYPE)
            {
                throw std::runtime_error("Only tetrahedral volume elements are accepted (element type " +
                                         std::to_string(SU2_TET_TYPE) + "), element type parsed = " + std::to_string(element_type) + "\n");
            }
            tet_connect.at(i) = t;
        }

        // Reading nodes
        ist >> tmp_string >> N_NODES;
        check_string_correctness(tmp_string, "NPOIN=");
        nodes.resize(N_NODES);

        for (Index i{0}; i < N_NODES; i++)
        {
            Vec3 node;
            ist >> node.x() >> node.y() >> node.z();
            nodes.at(i) = node;
        }

        // Reading boundary patches
        ist >> tmp_string >> N_PATCHES;
        check_string_correctness(tmp_string, "NMARK=");
        tri_patch_connect_list.reserve(N_PATCHES);

        for (Index i{0}; i < N_PATCHES; i++)
        {
            TriPatchConnect p;
            Index N_MARKER_ELEMENTS;
            ist >> tmp_string >> p.patch_name;
            check_string_correctness(tmp_string, "MARKER_TAG=");
            if (!config.input_file_contains_patch_name(p.patch_name))
                throw std::runtime_error("Patch with name '" + p.patch_name + "' is not named in input file\n");
            ist >> tmp_string >> N_MARKER_ELEMENTS;
            check_string_correctness(tmp_string, "MARKER_ELEMS=");
            p.triangles.resize(N_MARKER_ELEMENTS);
            for (Index j{0}; j < N_MARKER_ELEMENTS; j++)
            {
                TriConnect t;
                ist >> element_type >> t.a() >> t.b() >> t.c() >> element_num;
                if (element_type != SU2_TRI_TYPE)
                {
                    throw std::runtime_error("Only triangular surface elements are accepted (element type " +
                                             std::to_string(SU2_TRI_TYPE) + "), element type parsed = " + std::to_string(element_type) + "\n");
                }
                p.triangles.at(j) = t;
            }
            tri_patch_connect_list.push_back(p);
        }
    }

}