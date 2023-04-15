#include "../include/Grid.hpp"



void Grid::print_grid() const{

    cout << "CELLS:\n";
    for (Index i{0}; i< cells.size();i++) cout << i << ": " << cells.at(i) << endl;

    cout << "\n\nFACES:\n";
    for (Index i{0}; i< faces.size();i++) cout << i << ": " << faces.at(i) << endl;

}

void Grid::create_grid(Config& config){
    using namespace Geometry;


    //Create native mesh containers and read mesh file

    vector<Vec3> nodes; 
    vector<TetConnectivity> tet_connect; 
    vector<TriPatchConnectivity> tri_patch_connect;
    read_mesh(config.get_mesh_filename(), nodes, tet_connect, tri_patch_connect);

    patches.resize(tri_patch_connect.size());

    //The following procedure creates internal cells and faces

    Index N_INTERIOR_CELLS = tet_connect.size();
    cells.resize(N_INTERIOR_CELLS);
    Index N_NODES = nodes.size();

    TriConnectivity face_bound;
    vector<Triangle> face_triangles; //Used to assign properties to faces later
    Index next_ghost_index = N_NODES;
    
    for (Index i=0; i<N_INTERIOR_CELLS; i++){

        TetConnectivity tc_i = tet_connect.at(i);
        Tetrahedron tet{nodes.at(tc_i.a),nodes.at(tc_i.b), nodes.at(tc_i.c), nodes.at(tc_i.d)};

        //Adding a new Cell
        cells.emplace_back(tet.calc_volume(), tet.calc_centroid()); 

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

                face_bound = add_face_to_patches(face_ij, next_ghost_index, tri_patch_connect);
                Triangle tri_bound{nodes.at(face_bound.a), nodes.at(face_bound.b), nodes.at(face_bound.c)};                
                face_triangles.push_back(tri_bound);

                faces.emplace_back(i, next_ghost_index);

                next_ghost_index++;
            }
        }

        Index N_FACES = faces.size();
        Index N_TOTAL_CELLS = N_INTERIOR_CELLS + next_ghost_index;

        config.set_N_INTERIOR_CELLS(N_INTERIOR_CELLS);
        config.set_N_TOTAL_CELLS(N_TOTAL_CELLS);
        config.set_N_FACES(N_FACES);

        Index max_i{0}, max_j{0}; //For consistency checking

        //loop over all faces to assign area vector and centroid vectors
        assert(face_triangles.size() == N_FACES);
        for (Index ij{0}; ij<N_FACES; ij++){
            Face& face = faces.at(ij);
            Index i = face.i;
            Index j = face.j; 
            max_i = std::max(max_i, i);
            max_j = std::max(max_j, j);

            assign_face_properties(face, face_triangles.at(ij), cells.at(i).centroid, cells.at(j).centroid);

        }
    
        assert(max_i == N_INTERIOR_CELLS); //i should always belong to the interior cells
        assert(max_j == N_TOTAL_CELLS);
    }

    
}

void Grid::read_mesh(string mesh_filename, 
                   vector<Vec3>& nodes, 
                   vector<Geometry::TetConnectivity>& tet_connect, 
                   vector<Geometry::TriPatchConnectivity>& tri_patch_connect) const{
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
    string BC_type;
    Index N_surface_elements;    
    while (ist >> BC_type >> N_surface_elements){
        BoundaryType boundary_type = boundary_type_from_string.at(BC_type);
        
        TriConnectivity tri;
        TriPatchConnectivity p;
        p.boundary_type = boundary_type; 
        p.triangles.resize(N_surface_elements);
        for (Index i{0}; i<N_surface_elements; i++){
            ist >> tri.a >> tri.b >> tri.c;
            p.triangles.at(i) = tri;
        }
        tri_patch_connect.emplace_back(p);
    }
}


std::pair<Index,bool> Grid::find_neigbouring_cell(Index i, 
                                                Geometry::TriConnectivity face_ij, 
                                                const vector<Geometry::TetConnectivity>& tet_connect) const{
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
        if (face.i==i || face.i==j){
            assert(face.j==i || face.j==j);
            return true;        
        }
    return false;
}


Geometry::TriConnectivity Grid::add_face_to_patches(Geometry::TriConnectivity t_ij, 
                            Index ij, 
                            const vector<Geometry::TriPatchConnectivity>& tri_patch_connect){
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
                patches.at(patch_number).boundary_type = tri_patch_connect.at(patch_number).boundary_type;
                return t_patch;
            }
            
        }
    }
    std::cerr << "Error: the face triangle was not found among the boundary pathces\n";
    exit(1);             
}