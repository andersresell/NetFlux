#include "../include/Grid.hpp"

using namespace geom;

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
        cout << "patch type: " << (int)patch.boundary_type << 
        "\nFIRST FACE: " << patch.FIRST_FACE << 
        "\nN_FACES: " << patch.N_FACES << "\n\n";
    }

    cout << "\n\nFACE INDICES FROM CELLS:\n";
    for (Index i{0}; i < face_indices_from_cell.size(); i++)
    {
        cout << i << ": ";
        for (const auto &ind : face_indices_from_cell.at(i))
            cout << ind << " ";
        cout << endl;
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

void Grid::create_grid(Config &config)
{
    read_mesh(config.get_mesh_filename());
    Index N_NODES = nodes.size();
    
    create_interior();
    Index N_INTERIOR_CELLS = cells.size();
    Index N_INTERIOR_FACES = faces.size();

    create_boundaries(config);
    Index N_TOTAL_CELLS = cells.size();
    Index N_TOTAL_FACES = faces.size();

    config.set_grid_metrics(N_NODES, N_INTERIOR_CELLS, N_TOTAL_CELLS, N_INTERIOR_FACES, N_TOTAL_FACES);

    assign_geometry_properties(config);

    reorder_faces(config);

    shrink_vectors();
    face_triangles.clear(); //face_triangles are not needed anymore
}

void Grid::read_mesh(string mesh_filename)
{
    string extension = mesh_filename.substr(mesh_filename.find_last_of(".") + 1);

    if (extension == "c3d")
        read_c3d_mesh(mesh_filename);
    else if (extension == "su2")
        read_su2_mesh(mesh_filename);
    else
        FAIL_MSG("Error: Illegal mesh format");
}

void Grid::read_c3d_mesh(string mesh_filename)
{

    std::ifstream ist{mesh_filename};
    std::stringstream ss;

    FAIL_IF_MSG(!ist, "Error: couldn't open the mesh file: " + mesh_filename);

    string line, tmp;
    Index i, N_NODES, N_ELEMENTS;

    // READ N_NODES
    ist >> tmp >> N_NODES;
    FAIL_IF(tmp != "N_NODES");

    // READ N_ELEMENTS
    ist >> tmp;
    ist >> N_ELEMENTS;
    FAIL_IF(tmp != "N_ELEMENTS");

    nodes.resize(N_NODES);
    while (getline(ist, line))
    {
        if (line.size() != 0)
        {
            FAIL_IF(line != "#nodes");
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
            FAIL_IF(line != "#tetrahedra");
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
            FAIL_IF(line != "#patches");
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
        p.triangles.resize(N_surface_elements);
        for (Index i{0}; i < N_surface_elements; i++)
        {
            ist >> tri.a() >> tri.b() >> tri.c();
            p.triangles.at(i) = tri;
        }
        tri_patch_connect_list.emplace_back(p);
    }
}

void Grid::read_su2_mesh(string mesh_filename)
{

    std::ifstream ist{mesh_filename};
    std::stringstream ss;

    FAIL_IF_MSG(!ist, "Error: couldn't open the mesh file: " + mesh_filename);

    const ShortIndex SU2_TET_TYPE = 10, SU2_TRI_TYPE = 5;
    string tmp_string;
    ShortIndex tmp_int, element_type;
    Index N_NODES, N_ELEMENTS, N_PATCHES, element_num, first_element_num;

    ist >> tmp_string >> tmp_int;
    FAIL_IF(tmp_string != "NDIME=");
    FAIL_IF(tmp_int != 3);

    // Reading element connectivity. Only permitting tetrahedral type
    ist >> tmp_string >> N_ELEMENTS;
    FAIL_IF(tmp_string != "NELEM=");
    tet_connect.resize(N_ELEMENTS);

    for (Index i{0}; i < N_ELEMENTS; i++)
    {
        TetConnect t;
        ist >> element_type >> t.a() >> t.b() >> t.c() >> t.d() >> element_num;
        if (i == 0)
            first_element_num = element_num;
        FAIL_IF(element_type != SU2_TET_TYPE);
        if (i == N_ELEMENTS - 1)
            assert(element_num - first_element_num == N_ELEMENTS - 1);
        tet_connect.at(i) = t;
    }

    // Reading nodes
    ist >> tmp_string >> N_NODES;
    FAIL_IF(tmp_string != "NPOIN=");
    nodes.resize(N_NODES);

    for (Index i{0}; i < N_NODES; i++)
    {
        Vec3 node;
        ist >> node.x() >> node.y() >> node.z();
        nodes.at(i) = node;
    }

    // Reading boundary patches
    ist >> tmp_string >> N_PATCHES;
    FAIL_IF(tmp_string != "NMARK=");
    tri_patch_connect_list.reserve(N_PATCHES);

    for (Index i{0}; i < N_PATCHES; i++)
    {
        TriPatchConnect p;
        Index N_MARKER_ELEMENTS;
        ist >> tmp_string >> p.patch_name;
        FAIL_IF(tmp_string != "MARKER_TAG=");
        ist >> tmp_string >> N_MARKER_ELEMENTS;
        FAIL_IF(tmp_string != "MARKER_ELEMS=");
        p.triangles.resize(N_MARKER_ELEMENTS);
        for (Index j{0}; j < N_MARKER_ELEMENTS; j++)
        {
            TriConnect t;
            ist >> element_type >> t.a() >> t.b() >> t.c() >> element_num;
            FAIL_IF(element_type != SU2_TRI_TYPE);
            if (j == 0)
                first_element_num = element_num;
            if (j == N_MARKER_ELEMENTS - 1)
                FAIL_IF(element_num - first_element_num != N_MARKER_ELEMENTS - 1);
            p.triangles.at(j) = t;
        }
        tri_patch_connect_list.push_back(p);
    }
}

void Grid::create_interior(){
    
    Index N_TETS = tet_connect.size();
    Index N_INTERIOR_CELLS = N_TETS;
    
    cells.reserve(N_INTERIOR_CELLS + find_N_GHOST_cells());
    
    face_indices_from_cell.resize(N_INTERIOR_CELLS);
    Index N_NODES = nodes.size();

    for (Index i = 0; i < N_INTERIOR_CELLS; i++){

        TetConnect tc_i = tet_connect.at(i);
        Tetrahedron tet = tet_from_connect(tc_i);

        // Adding a new Cell
        cells.emplace_back(tet.calc_volume(), tet.calc_centroid());

        for (ShortIndex k{0}; k < N_TET_FACES; k++){

            TriConnect face_ij = tet_face_connectivity(tet_connect.at(i), k);
            Triangle tri = tri_from_connect(face_ij);

            std::pair<Index, bool> pair = find_neigbouring_cell(i, face_ij, tet_connect);
            bool neigbouring_cell_found = pair.second;

            if (neigbouring_cell_found)
            {
                Index j = pair.first;
                if (!face_ij_created(i, j))
                {
                    // Create new face
                    faces.emplace_back(i, j);
                    add_face_to_cell_i(i, j);
                    face_triangles.push_back(tri);
                }
            }
        }
    }
}

void Grid::create_boundaries(const Config& config){

    const Index N_PATCHES = tri_patch_connect_list.size(); 
    
    faces.reserve(faces.size() + find_N_GHOST_cells());
    
    for (const auto &tpc : tri_patch_connect_list){
 
        Patch p;
        p.boundary_type = config.get_boundary_type(tpc.patch_name);
        p.FIRST_FACE = faces.size();
        p.N_FACES = tpc.triangles.size(); 
        patches.push_back(p);

        for (const TriConnect& tc : tpc.triangles){
            Index i = find_boundary_face_owner(tc);
            Index j_ghost = cells.size();
            cells.emplace_back(); //Adding ghost cell
            faces.emplace_back(i, j_ghost);
            face_triangles.push_back(tri_from_connect(tc));
        }   
    }
}


void Grid::assign_geometry_properties(const Config& config){

    Index max_i{0}, max_j{0}; // For consistency checking

    // loop over all faces to assign area vector and centroid vectors
    assert(face_triangles.size() == config.get_N_TOTAL_FACES());

    Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();
    Index N_TOTAL_CELLS = config.get_N_TOTAL_CELLS();

    for (Index ij{0}; ij < config.get_N_TOTAL_FACES(); ij++)
    {
        Face &face = faces.at(ij);
        Index i = face.i;
        Index j = face.j;
        max_i = std::max(max_i, i);
        max_j = std::max(max_j, j);


        assert(i < N_INTERIOR_CELLS); // i should never belong to a ghost cell (normal pointing outwards)

        if (j >= N_INTERIOR_CELLS)
        {
            cells.at(j).centroid = calc_ghost_centroid(cells.at(i).centroid, face_triangles.at(ij));
            cells.at(j).cell_volume = 0;
        }

        assign_face_properties(face, face_triangles.at(ij), cells.at(i).centroid, cells.at(j).centroid);
    }

    assert(max_i == N_INTERIOR_CELLS - 1); // i should always belong to the interior cells
    assert(max_j == N_TOTAL_CELLS - 1);

    assert(N_TOTAL_CELLS == cells.size());
}

Index Grid::find_boundary_face_owner(TriConnect tc){
    /*Find the interior cell of the boundary face with a given connectivity.
    It works by comparing the tetrahedral element connectivity to tc. If all {a b c} in tc is contained in {a b c d} 
    of the element, then element i is the correct owner*/
    
    TriConnect face_ij;
    TetConnect tet_i;
    array<Index, N_TRI_NODES> boundary_compare, ij_compare;
    bool all_nodes_equal;

    boundary_compare = {tc.a(), tc.b(), tc.c()};
    std::sort(boundary_compare.begin(), boundary_compare.end());

    for (Index i{0}; i<tet_connect.size(); i++){
        tet_i = tet_connect.at(i);

        for (ShortIndex k{0}; k < N_TET_FACES; k++){
            
            face_ij = tet_face_connectivity(tet_i, k);
            ij_compare = {face_ij.a(), face_ij.b(), face_ij.c()};
            std::sort(ij_compare.begin(), ij_compare.end());

            all_nodes_equal = true;
            for (ShortIndex l{0}; l < N_TRI_NODES; l++)
                if (ij_compare.at(l) != boundary_compare.at(l)){
                    all_nodes_equal = false;
                    break;
                }

            if (all_nodes_equal)
                return i;
        }
    }
    FAIL_MSG("Error, no matching cell was found for the provided boundary face\n");
}


void Grid::reorder_faces(const Config& config){
    /*Sorts the faces based on the indices of the neighbour cells. We want them to be sorted in increasing order,
    with the owner cell i being prioritized first and the neighbour cell j second. As an example,
    the ordering {(0,1), (0,3), (0,2), (1,1)} should be changed to {(0,1), (0,2), (0,3), (1,1)}
    The way the grid is constructed this is allready achieved for the owner cell i. However j might need sorting. 
    The comparison operator < for Face is defined for this purpose.
    The interior and each boundary patch is sorted separately.*/
    Index N_INTERIOR_FACES = config.get_N_INTERIOR_FACES();
    std::sort(faces.begin(), faces.begin() + N_INTERIOR_FACES);

    for (Patch& patch : patches){
        std::sort(faces.begin() + patch.FIRST_FACE, faces.begin() + patch.FIRST_FACE + patch.N_FACES);
    }
}

std::pair<Index, bool> Grid::find_neigbouring_cell(Index i,
                                                   TriConnect face_ij,
                                                   const Vector<TetConnect> &tet_connect) const
{
    /*Finds the neighbouring cell j of face ij of cell i. If no neigbour exist (boundary), it returns false.
    The function works by looping over all other cells. For each cell it loops over all faces.
    It then compares the face nodes to the face nodes of face ij (both sorted). If they are equal the cell is found*/

    Index N_CELLS = tet_connect.size();
    array<Index, N_TRI_NODES> i_compare, j_compare;
    TriConnect face_ji;
    TetConnect tet_j;
    bool all_nodes_equal;

    i_compare = {face_ij.a(), face_ij.b(), face_ij.c()};
    std::sort(i_compare.begin(), i_compare.end());

    for (Index j{0}; j < N_CELLS; j++)
    {
        if (i == j)
            continue;

        tet_j = tet_connect.at(j);

        for (ShortIndex k{0}; k < N_TET_FACES; k++)
        {
            face_ji = tet_face_connectivity(tet_j, k);
            j_compare = {face_ji.a(), face_ji.b(), face_ji.c()};
            std::sort(j_compare.begin(), j_compare.end());

            all_nodes_equal = true;
            for (ShortIndex l{0}; l < N_TRI_NODES; l++)
                if (i_compare.at(l) != j_compare.at(l))
                    all_nodes_equal = false;

            if (all_nodes_equal)
                return {j, true};
        }
    }
    return {N_CELLS, false};
}

bool Grid::face_ij_created(Index i, Index j) const{
    for (const Face &face : faces)
        if ((face.i == i || face.i == j) && (face.j == i || face.j == j))
            return true;
    
    return false;
}

Tetrahedron Grid::tet_from_connect(const TetConnect &tc) const{
    return Tetrahedron(nodes.at(tc.a()), nodes.at(tc.b()), nodes.at(tc.c()), nodes.at(tc.d()));
}

Triangle Grid::tri_from_connect(const TriConnect &tc) const{
    return Triangle(nodes.at(tc.a()), nodes.at(tc.b()), nodes.at(tc.c()));
}

void Grid::add_face_to_cell_i(Index i, Index ij){
    face_indices_from_cell.at(i).push_back(ij);
}


Index Grid::find_N_GHOST_cells(){
    Index N_GHOST{0};
    for (const auto& tpc : tri_patch_connect_list)
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
    cells.shrink_to_fit();
    faces.shrink_to_fit();
    patches.shrink_to_fit();

    for (auto &indices : face_indices_from_cell)
        indices.shrink_to_fit();
    face_indices_from_cell.shrink_to_fit();
}