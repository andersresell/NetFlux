#include "../include/Grid.hpp"



void Grid::print_grid(const Config& config) const{

    cout << "CELLS:\n";
    for (Index i{0}; i< cells.size();i++) {
        if (i >= config.get_N_INTERIOR_CELLS())
            cout << "GHOST, ";
        cout << i << ": " << cells.at(i) << endl;
    }
    
    cout << "\n\nFACES:\n";
    for (Index i{0}; i< faces.size();i++) cout << i << ": " << faces.at(i) << endl;
    
    cout << "\n\nPATCHES:\n";
    for (const auto& patch : patches){
        cout << "patch type: " << (int)patch.boundary_type << "\nFace indices:\n";
        for (Index ij : patch.boundary_face_indices) cout << ij << " ";
        cout << "\n\n";
    }
}

void Grid::print_native_mesh() const{
    using namespace Geometry;
    cout << "NODES:\n";
    for (Index i{0}; i<nodes.size();i++) {
        cout << i << ": " << horizontal_string_Vec3(nodes.at(i)) << endl;
    }
    cout << "\n\nTET CONNECTIVITY:\n";
    for (Index i{0}; i<tet_connect.size();i++) {
        TetConnectivity t = tet_connect.at(i);
        cout << i << ": " << t.a << ", " << t.b << ", " << t.c << ", " << t.d << endl;
    }
    cout << "\n\nTRIANGLE PATCH CONNECTIVITY:\n";
    for (const auto& tpc : tri_patch_connect){
        cout << "BC type: " << tpc.patch_name << endl;
        for (Index i{0}; i<tpc.triangles.size();i++) {
            TriConnectivity t = tpc.triangles.at(i);
            cout << i << ": " << t.a << ", " << t.b << ", " << t.c << endl;
        }
        cout << "\n\n";
    }
    
}

void Grid::create_grid(Config& config){
    using namespace Geometry;

    /*Read mesh file. This populates the:
    - nodes -> native mesh nodes
    - tet_connect -> tetrahedral elements connectivity
    - tri_patch_connect -> patches of boundary triangles and respective BC types
     */
    read_mesh(config.get_mesh_filename());
    

    //Copying boundary condition from tri_patch_connect to the patches
    assign_patch_BC(config);
                
    print_native_mesh();

    Index N_TETS = tet_connect.size();
    Index N_INTERIOR_CELLS = N_TETS;
    cells.resize(N_INTERIOR_CELLS);
    Index N_NODES = nodes.size();

    TriConnectivity face_bound;
    vector<Triangle> face_triangles; //Used to assign properties to faces later
    Index next_ghost_index = N_INTERIOR_CELLS;
    


    for (Index i=0; i<N_INTERIOR_CELLS; i++){

        TetConnectivity tc_i = tet_connect.at(i);
        Tetrahedron tet{nodes.at(tc_i.a),nodes.at(tc_i.b), nodes.at(tc_i.c), nodes.at(tc_i.d)};

        //Adding a new Cell
        cells.at(i) = Cell{tet.calc_volume(), tet.calc_centroid()};

        for (ShortIndex k{0}; k<N_TET_FACES; k++){

            TriConnectivity face_ij = tet_face_connectivity(tet_connect.at(i), k); 
            Triangle tri{nodes.at(face_ij.a), nodes.at(face_ij.b), nodes.at(face_ij.c)};
            std::pair<Index,bool> pair = find_neigbouring_cell(i, face_ij, tet_connect);
            bool neigbouring_cell_found = pair.second;

            if (neigbouring_cell_found){
                Index j = pair.first;
                if (!face_ij_created(i, j)){
                    //Create new face
                    faces.emplace_back(i, j);
                    face_triangles.push_back(tri);
                }
            } else{ //No neigbour found. Marking a boundary face
                Index next_ij = faces.size();
                face_bound = add_face_to_patches(face_ij, next_ij, tri_patch_connect);
                Triangle tri_bound{nodes.at(face_bound.a), nodes.at(face_bound.b), nodes.at(face_bound.c)};                
                face_triangles.push_back(tri_bound);

                faces.emplace_back(i, next_ghost_index);
                cells.resize(next_ghost_index+1);
                cells.at(next_ghost_index) =  Cell{};
                next_ghost_index++;
            }
        }
    }

    Index N_FACES = faces.size();
    Index N_TOTAL_CELLS = next_ghost_index;

    config.set_grid_data(N_NODES,N_INTERIOR_CELLS, N_TOTAL_CELLS, N_FACES);

    Index max_i{0}, max_j{0}; //For consistency checking

    //loop over all faces to assign area vector and centroid vectors
    assert(face_triangles.size() == N_FACES);
    for (Index ij{0}; ij<N_FACES; ij++){
        Face& face = faces.at(ij);
        Index i = face.i;
        Index j = face.j; 
        max_i = std::max(max_i, i);
        max_j = std::max(max_j, j);

        assert(i < N_INTERIOR_CELLS); //i should never belong to a ghost cell (normal pointing outwards)
        //if j is a ghost cell its centoid is calculated fromt he interior 
        //Vec3 centroid_j = (j >= N_INTERIOR_CELLS) ? calc_ghost_centroid(cells.at(i).centroid, face_triangles.at(ij)) : cells.at(j).centroid;
        
        //assign_face_properties(face, face_triangles.at(ij), cells.at(i).centroid, centroid_j);
        
        if (j >= N_INTERIOR_CELLS){
            //Vec3 centroid_j = calc_ghost_centroid(cells.at(i).centroid, face_triangles.at(ij));
            cells.at(j).centroid = calc_ghost_centroid(cells.at(i).centroid, face_triangles.at(ij));
            cells.at(j).cell_volume = 0;       
        }
        
        assign_face_properties(face, face_triangles.at(ij), cells.at(i).centroid, cells.at(j).centroid);

    }

    assert(max_i == N_INTERIOR_CELLS - 1); //i should always belong to the interior cells
    assert(max_j == N_TOTAL_CELLS - 1);

    assert(N_TOTAL_CELLS == cells.size());
    
}

void Grid::read_mesh(string mesh_filename){
    string extension = mesh_filename.substr(mesh_filename.find_last_of(".")+1);
    
    if (extension == "c3d") read_c3d_mesh(mesh_filename);
    else if (extension == "su2") read_su2_mesh(mesh_filename);
    else {
        std::cerr << "Error: Illegal mesh format\n";
        exit(1);
    }
}

void Grid::read_c3d_mesh(string mesh_filename){
    using namespace Geometry;

    std::ifstream ist{mesh_filename};
    std::stringstream ss;

    if (!ist) {
        std::cerr << "Error: couldn't open the mesh file: "+mesh_filename;
        exit(1);
    }
    
    string line, tmp;
    Index i, N_NODES, N_ELEMENTS;

    //READ N_NODES
    ist >> tmp >> N_NODES;
    assert(tmp == "N_NODES");
    
    //READ N_ELEMENTS
    ist >> tmp; ist >> N_ELEMENTS;
    assert(tmp == "N_ELEMENTS");

    nodes.resize(N_NODES);
    while (getline(ist, line)) {
        if (line.size() != 0){
            assert(line == "#nodes");
            break;
        }
    }
    //reading nodes
    for (Index i{0}; i<N_NODES; i++){
        Vec3 node;
        ist >> node.x() >> node.y() >> node.z();
        nodes.at(i) = node;
    }
    
    tet_connect.resize(N_ELEMENTS);
    while (getline(ist, line)) {
        if (line.size() != 0){
            assert(line == "#tetrahedra");
            break;
        } 
    }
    //reading elements (assuming only tetrahedra for now)
    for (Index i{0}; i<N_ELEMENTS; i++){
        TetConnectivity t;
        ist >> t.a >> t.b >> t.c >> t.d;
        tet_connect.at(i) = t;
    }

    while (getline(ist, line)) {
        if (line.size() != 0){
            assert(line == "#patches");
            break;
        } 
    }

    //Reading boundary patches
    string patch_name;
    Index N_surface_elements;    
    while (ist >> patch_name >> N_surface_elements){
        TriConnectivity tri;
        TriPatchConnectivity p;
        p.patch_name = patch_name;
        p.triangles.resize(N_surface_elements);
        for (Index i{0}; i<N_surface_elements; i++){
            ist >> tri.a >> tri.b >> tri.c;
            p.triangles.at(i) = tri;
        }
        tri_patch_connect.emplace_back(p);
    }
}

void Grid::read_su2_mesh(string mesh_filename){
    using namespace Geometry;

    std::ifstream ist{mesh_filename};
    std::stringstream ss;

    if (!ist) {
        std::cerr << "Error: couldn't open the mesh file: "+mesh_filename;
        exit(1);
    }
    const ShortIndex SU2_TET_TYPE = 10, SU2_TRI_TYPE = 5;
    string tmp_string;
    ShortIndex tmp_int, element_type; 
    Index N_NODES, N_ELEMENTS, N_PATCHES, element_num, first_element_num;

    ist >> tmp_string >> tmp_int;
    assert(tmp_string == "NDIME=");
    assert(tmp_int == 3);

    //Reading element connectivity. Only permitting tetrahedral type
    ist >> tmp_string >> N_ELEMENTS;
    assert(tmp_string == "NELEM=");
    tet_connect.resize(N_ELEMENTS);
    
    for (Index i{0}; i<N_ELEMENTS; i++){ 
        TetConnectivity t;
        ist >> element_type >> t.a >> t.b >> t.c >> t.d >> element_num;
        if (i == 0) first_element_num = element_num;    
        assert(element_type == SU2_TET_TYPE);
        if (i == N_ELEMENTS-1) assert(element_num - first_element_num == N_ELEMENTS-1);
        tet_connect.at(i) = t;
    }

    //Reading nodes
    ist >> tmp_string >> N_NODES;
    assert(tmp_string == "NPOIN=");
    nodes.resize(N_NODES);

    for (Index i{0}; i<N_NODES; i++){
        Vec3 node;
        ist >> node.x() >> node.y() >> node.z();
        nodes.at(i) = node;
    }

    //Reading boundary patches
    ist >> tmp_string >> N_PATCHES;
    assert(tmp_string == "NMARK=");
    tri_patch_connect.reserve(N_PATCHES);

    for (Index i{0}; i<N_PATCHES; i++){
        TriPatchConnectivity p;
        Index N_MARKER_ELEMENTS;
        ist >> tmp_string >> p.patch_name;
        assert(tmp_string == "MARKER_TAG=");
        ist >> tmp_string >> N_MARKER_ELEMENTS;
        assert(tmp_string == "MARKER_ELEMS=");
        p.triangles.resize(N_MARKER_ELEMENTS);
        for (Index j{0}; j<N_MARKER_ELEMENTS; j++){
            TriConnectivity t;   
            ist >> element_type >> t.a  >> t.b >> t.c >> element_num;
            assert(element_type == SU2_TRI_TYPE);
            if (j == 0) first_element_num = element_num; 
            if (j == N_MARKER_ELEMENTS-1) assert(element_num - first_element_num == N_MARKER_ELEMENTS-1);
            p.triangles.at(j) = t;
        }
        tri_patch_connect.push_back(p);
    }

}

 
std::pair<Index,bool> Grid::find_neigbouring_cell(Index i, 
                                                Geometry::TriConnectivity face_ij, 
                                                const vector<Geometry::TetConnectivity>& tet_connect) const{
    /*Finds the neighbouring cell j of face ij of cell i. If no neigbour exist (boundary), it returns false.
    The function works by looping over all other cells. For each cell it loops over all faces. 
    It then compares the face nodes to the face nodes of face ij (both sorted). If they are equal the cell is found*/
    using namespace Geometry;
    Index N_CELLS = tet_connect.size();
    array<Index, N_TRI_NODES> i_compare, j_compare;
    TriConnectivity face_ji;
    TetConnectivity tet_j;
    bool all_nodes_equal;

    i_compare = {face_ij.a, face_ij.b, face_ij.c};
    std::sort(i_compare.begin(), i_compare.end());

    for (Index j{0}; j<N_CELLS; j++){
        if (i==j) continue;

        tet_j = tet_connect.at(j);

        for (ShortIndex k{0}; k<N_TET_FACES; k++){
            face_ji = tet_face_connectivity(tet_j, k);
            j_compare = {face_ji.a, face_ji.b, face_ji.c};
            std::sort(j_compare.begin(), j_compare.end());
            
            all_nodes_equal = true;
            for (ShortIndex l{0}; l<N_TRI_NODES; l++)
                if (i_compare.at(l) != j_compare.at(l)) all_nodes_equal = false;
            
            if (all_nodes_equal) return {j, true};
        }

    }
    return {N_CELLS, false};
}


bool Grid::face_ij_created(Index i, Index j) const{
    for (const Face& face : faces)
        if ((face.i==i || face.i==j) && (face.j==i || face.j==j)){
            return true;        
        }
    return false;
}


Geometry::TriConnectivity Grid::add_face_to_patches(Geometry::TriConnectivity t_ij, 
                            Index ij, 
                            const vector<Geometry::TriPatchConnectivity>& tri_patch_connect){
    /*Function adds the index of face ij to the patch where the triangle t_ij belongs*/
    using namespace Geometry;
    
    array<Index,N_TRI_NODES> ij_compare{t_ij.a, t_ij.b, t_ij.c}, patch_compare;
    std::sort(ij_compare.begin(), ij_compare.end());
    
    assert(tri_patch_connect.size() == patches.size());

    for (Index patch_number{0}; patch_number<patches.size(); patch_number++){
        for (const TriConnectivity& t_patch : tri_patch_connect.at(patch_number).triangles){
            
            patch_compare = {t_patch.a, t_patch.b, t_patch.c};
            std::sort(patch_compare.begin(), patch_compare.end());
            
            if (arrays_equal(ij_compare, patch_compare))
            {
                //Correct triangle found
                patches.at(patch_number).boundary_face_indices.push_back(ij);
                return t_patch;
            }
            
        }
    }
    std::cerr << "Error: the face triangle was not found among the boundary pathces.\n" 
        << "Some boundary triangles might not have been assigned patches\n";
    exit(1);             
}

 void Grid::assign_patch_BC(const Config& config){
    patches.resize(tri_patch_connect.size());
    for (Index i{0}; i<tri_patch_connect.size(); i++){
        string patch_name = tri_patch_connect.at(i).patch_name;
        patches.at(i).boundary_type = config.get_boundary_type(patch_name);
    }
 } 
    