#include "../include/Grid.hpp"

using namespace geom;

Grid::Grid(Config &config)
{
    read_mesh(config);

    try
    {
        create_grid(config);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Error creating grid:\n" + string(e.what()));
    }
}

void Grid::read_mesh(const Config &config)
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

void Grid::read_netflux_mesh(const Config &config)
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

void Grid::read_su2_mesh(const Config &config)
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

void Grid::create_grid(Config &config)
{
    cout << "Creating grid from native mesh\n";

    /*--------------------------------------------------------------------
    Step 1: Create all cells from elements, this is straight forward
    --------------------------------------------------------------------*/

    cells.reserve(tet_connect.size() + find_N_GHOST_cells());
    for (Index i{0}; i < tet_connect.size(); i++)
    {
        TetConnect tc_i = tet_connect[i];
        Tetrahedron tet = tet_from_connect(tc_i);
        cells.cell_volumes.emplace_back(tet.calc_volume());
        cells.centroids.emplace_back(tet.calc_centroid());
    }

    /*--------------------------------------------------------------------
    Step 2: Associate the face nodes with its two neighboring cells.
    --------------------------------------------------------------------*/

    Vector<Index> face_indices;
    map<SortedTriConnect, pair<Index, long int>> face_to_cells;
    constexpr int CELL_NOT_YET_ASSIGNED{-1};

    for (Index cell_index{0}; cell_index < tet_connect.size(); cell_index++)
    {

        for (ShortIndex k{0}; k < N_TET_FACES; k++)
        {
            TriConnect tc{tet_face_connectivity(tet_connect[cell_index], k)};

            // Get the sorted node indices of face k
            SortedTriConnect face_ij{tc};

            assert(face_to_cells.count(face_ij) <= 1);
            if (face_to_cells.count(face_ij) == 0)
            {
                // Discovered a new face
                face_to_cells.insert({face_ij, {cell_index, CELL_NOT_YET_ASSIGNED}});
            }
            else
            {
                // Face allready discovered by previous cell
                face_to_cells.at(face_ij).second = cell_index;
            }
        }
    }
    Vector<Triangle> face_triangles;
    face_triangles.reserve(face_to_cells.size());
    faces.resize(face_to_cells.size());

    /*--------------------------------------------------------------------
    Step 3: Adding internal faces to the map. (Boundary faces are added
    later, this is to get the correct grouping of patches). Face triangles
    (used for calculating geometry properties) are created simultaneously
    as faces, to get the correct ordering.
    --------------------------------------------------------------------*/

    for (const auto &face : face_to_cells)
    {
        Index cell_i = face.second.first;
        long int cell_j = face.second.second;
        if (cell_j != CELL_NOT_YET_ASSIGNED)
        { // Only add internal faces
            assert(cell_i < cell_j);
            faces.cell_indices.emplace_back(cell_i, cell_j);
            face_triangles.emplace_back(tri_from_connect(face.first));
        }
    }
    assert(face_triangles.size() == faces.size());

    /*--------------------------------------------------------------------
    Step 4: Create the faces at the boundaries and ghost cells.
    Also add the face triangles at the boundaries.
    --------------------------------------------------------------------*/
    cout << "Creating ghost cells..\n";
    for (const auto &tpc : tri_patch_connect_list)
    {

        Patch p;
        p.boundary_type = config.get_boundary_type(tpc.patch_name);
        p.FIRST_FACE = faces.size();
        p.N_FACES = tpc.triangles.size();
        patches.push_back(p);

        for (const TriConnect &tc : tpc.triangles)
        {
            SortedTriConnect face_ij{tc};
            Index cell_j_ghost = cells.size();
            assert(face_to_cells.at(face_ij).second == CELL_NOT_YET_ASSIGNED);
            face_to_cells.at(face_ij).second = cell_j_ghost; // this won't be used, so strictly unnecessary
            Index cell_i_domain = face_to_cells.at(face_ij).first;
            assert(cell_i_domain < cell_j_ghost);
            assert(cell_j_ghost >= tet_connect.size());
            faces.cell_indices.emplace_back(cell_i_domain, cell_j_ghost);

            face_triangles.emplace_back(tri_from_connect(tc));
            cells.add_empty(); // Adding ghost cell
        }
    }
    assert(face_triangles.size() == faces.size());

    /*-------------------------------------------------------------------
        Set some grid metrics in the config object
    --------------------------------------------------------------------*/
    Index N_NODES = nodes.size();
    Index N_INTERIOR_CELLS = tet_connect.size();
    Index N_TOTAL_CELLS = cells.size();
    Index N_INTERIOR_FACES = faces.size() - find_N_GHOST_cells();
    Index N_TOTAL_FACES = faces.size();
    config.set_grid_metrics(N_NODES, N_INTERIOR_CELLS, N_TOTAL_CELLS, N_INTERIOR_FACES, N_TOTAL_FACES);

    /*--------------------------------------------------------------------
    Sort faces so that the interior faces appear first with the owner index
    i allways being less than neigbour index j. The same logic is applied
    patch-wise to the boundaries
    --------------------------------------------------------------------*/
    cout << "Reorder faces..\n";
    reorder_faces(config);

    /*--------------------------------------------------------------------
    Assigning geometrical properties to faces and ghost cells
    --------------------------------------------------------------------*/
    cout << "Assign geometrical properties..\n";
    assign_geometry_properties(config, face_triangles);

    /*--------------------------------------------------------------------
    Reducing allocated memory
    --------------------------------------------------------------------*/
    shrink_vectors();

    cout << "Computational grid has been created.\n";
}

void Grid::assign_geometry_properties(const Config &config, const Vector<Triangle> &face_triangles)
{

    Index max_i{0}, max_j{0}; // For consistency checking

    // loop over all faces to assign area vector and centroid vectors
    assert(face_triangles.size() == config.get_N_TOTAL_FACES());

    Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();

    for (Index ij{0}; ij < config.get_N_TOTAL_FACES(); ij++)
    {
        Index i = faces.get_cell_i(ij);
        Index j = faces.get_cell_j(ij);
        max_i = std::max(max_i, i);
        max_j = std::max(max_j, j);

        assert(i < N_INTERIOR_CELLS); // i should never belong to a ghost cell (normal pointing outwards)

        if (j >= N_INTERIOR_CELLS)
        {
            cells.cell_volumes[j] = 0.0;
            cells.centroids[j] = calc_ghost_centroid(cells.centroids[i], face_triangles[ij]);
        }

        assign_face_properties(faces.normal_areas[ij],
                               faces.centroid_to_face_i[ij],
                               faces.centroid_to_face_j[ij],
                               face_triangles[ij],
                               cells.centroids[i],
                               cells.centroids[j]);
    }

    assert(max_i == config.get_N_INTERIOR_CELLS() - 1);
    assert(max_j == config.get_N_TOTAL_CELLS() - 1);
    assert(config.get_N_TOTAL_CELLS() == cells.size());
}

// Index Grid::find_boundary_face_owner(TriConnect tc)
// {
//     /*Find the interior cell of the boundary face with a given connectivity.
//     It works by comparing the tetrahedral element connectivity to tc. If all {a b c} in tc is contained in {a b c d}
//     of the element, then element i is the correct owner*/

//     TriConnect face_ij;
//     TetConnect tet_i;
//     array<Index, N_TRI_NODES> boundary_compare, ij_compare;
//     bool all_nodes_equal;

//     boundary_compare = {tc.a(), tc.b(), tc.c()};
//     std::sort(boundary_compare.begin(), boundary_compare.end());

//     for (Index i{0}; i < tet_connect.size(); i++)
//     {
//         tet_i = tet_connect.at(i);

//         for (ShortIndex k{0}; k < N_TET_FACES; k++)
//         {

//             face_ij = tet_face_connectivity(tet_i, k);
//             ij_compare = {face_ij.a(), face_ij.b(), face_ij.c()};
//             std::sort(ij_compare.begin(), ij_compare.end());

//             all_nodes_equal = true;
//             for (ShortIndex l{0}; l < N_TRI_NODES; l++)
//                 if (ij_compare.at(l) != boundary_compare.at(l))
//                 {
//                     all_nodes_equal = false;
//                     break;
//                 }

//             if (all_nodes_equal)
//                 return i;
//         }
//     }
//     FAIL_MSG("Error, no matching cell was found for the provided boundary face\n");
// }

void Grid::reorder_faces(const Config &config)
{
    /*Sorts the faces based on the indices of the neighbour cells. We want them to be sorted in increasing order,
    with the owner cell i being prioritized first and the neighbour cell j second. As an example,
    the ordering {(0,1), (0,3), (0,2), (1,1)} should be changed to {(0,1), (0,2), (0,3), (1,1)}
    The way the grid is constructed this is allready achieved for the owner cell i. However j might need sorting.
    The comparison operator < for Face is defined for this purpose.
    The interior and each boundary patch is sorted separately.*/
    Index N_INTERIOR_FACES = config.get_N_INTERIOR_FACES();
    std::sort(faces.cell_indices.begin(), faces.cell_indices.begin() + N_INTERIOR_FACES);

    for (Patch &patch : patches)
    {
        std::sort(faces.cell_indices.begin() + patch.FIRST_FACE, faces.cell_indices.begin() + patch.FIRST_FACE + patch.N_FACES);
    }
}

// std::pair<Index, bool> Grid::find_neigbouring_cell(Index i,
//                                                    TriConnect face_ij,
//                                                    const Vector<TetConnect> &tet_connect) const
// {
//     /*Finds the neighbouring cell j of face ij of cell i. If no neigbour exist (boundary), it returns false.
//     The function works by looping over all other cells. For each cell it loops over all faces.
//     It then compares the face nodes to the face nodes of face ij (both sorted). If they are equal the cell is found*/

//     Index N_CELLS = tet_connect.size();
//     array<Index, N_TRI_NODES> i_compare, j_compare;
//     TriConnect face_ji;
//     TetConnect tet_j;
//     bool all_nodes_equal;

//     i_compare = {face_ij.a(), face_ij.b(), face_ij.c()};
//     std::sort(i_compare.begin(), i_compare.end());

//     for (Index j{0}; j < N_CELLS; j++)
//     {
//         if (i == j)
//             continue;

//         tet_j = tet_connect.at(j);

//         for (ShortIndex k{0}; k < N_TET_FACES; k++)
//         {
//             face_ji = tet_face_connectivity(tet_j, k);
//             j_compare = {face_ji.a(), face_ji.b(), face_ji.c()};
//             std::sort(j_compare.begin(), j_compare.end());

//             all_nodes_equal = true;
//             for (ShortIndex l{0}; l < N_TRI_NODES; l++)
//                 if (i_compare.at(l) != j_compare.at(l))
//                     all_nodes_equal = false;

//             if (all_nodes_equal)
//                 return {j, true};
//         }
//     }
//     return {N_CELLS, false};
// }

// bool Grid::face_ij_created(Index i, Index j) const
// {
//     for (const Face &face : faces)
//         if ((face.i == i || face.i == j) && (face.j == i || face.j == j))
//             return true;

//     return false;
// }

Tetrahedron Grid::tet_from_connect(const TetConnect &tc) const
{
    return Tetrahedron(nodes.at(tc.a()), nodes.at(tc.b()), nodes.at(tc.c()), nodes.at(tc.d()));
}

Triangle Grid::tri_from_connect(const TriConnect &tc) const
{
    return Triangle(nodes.at(tc.a()), nodes.at(tc.b()), nodes.at(tc.c()));
}

Index Grid::find_N_GHOST_cells()
{
    Index N_GHOST{0};
    for (const auto &tpc : tri_patch_connect_list)
        N_GHOST += tpc.triangles.size();
    return N_GHOST;
}

void Grid::shrink_vectors()
{
    nodes.shrink_to_fit();
    tet_connect.shrink_to_fit();

    for (auto &tpc : tri_patch_connect_list)
        tpc.triangles.shrink_to_fit();

    tri_patch_connect_list.shrink_to_fit();
    cells.cell_volumes.shrink_to_fit();
    cells.centroids.shrink_to_fit();
    faces.cell_indices.shrink_to_fit();
    faces.centroid_to_face_i.shrink_to_fit();
    faces.centroid_to_face_j.shrink_to_fit();
    patches.shrink_to_fit();
}

void Grid::print_grid(const Config &config) const
{

    cout << "CELLS:\n";
    for (Index i{0}; i < cells.size(); i++)
    {
        if (i >= config.get_N_INTERIOR_CELLS())
            cout << "GHOST, ";
        cout << i << ": " << cells.at(i) << endl;
    }

    cout << "\n\nFACES:\n";
    for (Index i{0}; i < faces.size(); i++)
        cout << i << ": " << faces.at(i) << endl;

    cout << "\n\nPATCHES:\n";
    for (const auto &patch : patches)
    {
        cout << "patch type: " << (int)patch.boundary_type << "\nFIRST FACE: " << patch.FIRST_FACE << "\nN_FACES: " << patch.N_FACES << "\n\n";
    }
}

void Grid::print_native_mesh() const
{
    cout << "NODES:\n";
    for (Index i{0}; i < nodes.size(); i++)
    {
        cout << i << ": " << horizontal_string_Vec3(nodes.at(i)) << endl;
    }
    cout << "\n\nTET CONNECTIVITY:\n";
    for (Index i{0}; i < tet_connect.size(); i++)
    {
        TetConnect t = tet_connect.at(i);
        cout << i << ": " << t << endl;
    }
    cout << "\n\nTRIANGLE PATCH CONNECTIVITY:\n";
    for (const auto &tpc : tri_patch_connect_list)
    {
        cout << "BC type: " << tpc.patch_name << endl;
        for (Index i{0}; i < tpc.triangles.size(); i++)
        {
            TriConnect t = tpc.triangles.at(i);
            cout << i << ": " << t << endl;
        }
        cout << "\n\n";
    }
}