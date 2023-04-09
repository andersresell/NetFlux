#include "../include/Grid.hpp"

Grid::Grid(const Config& config){

}


void Grid::create_grid(const Config& config){

    //Create native mesh containers and read mesh file

    vector<Vec3> nodes; 
    vector<Geometry::TetConnectivity> tet_connect; 
    vector<Geometry::TriPatchConnectivity> tri_patch_connect;
    read_mesh(config.get_mesh_filename(), nodes, tet_connect, tri_patch_connect);

    //Then fill containers
    Index N_CELLS = tet_connect.size();
    cells.resize(N_CELLS);
    Index N_NODES = nodes.size();
    
    array<Index, 4> cell_i_compare, cell_j_compare;


    for (Index cell_i=0; cell_i<N_CELLS; cell_i++){
        Geometry::TetConnectivity tc_i = tet_connect.at(cell_i);
        Vec3 a = nodes.at(tc_i.a);
        Vec3 b = nodes.at(tc_i.b);
        Vec3 c = nodes.at(tc_i.c);
        Vec3 d = nodes.at(tc_i.d);
        Geometry::Tetrahedron tet{a,b,c,d};

        //Adding a new Cell
        cells.emplace_back(tet.calc_volume(), tet.calc_centroid()); 

        cell_i_compare = {tc_i.a, tc_i.b, tc_i.c, tc_i.d};
        std::sort(cell_i_compare.begin(), cell_i_compare.end());
        std::set<Index> set_i(cell_i_compare.begin(), cell_i_compare.end());
        bool neigbour_found = false;
        for (Index cell_j=0; cell_j<N_CELLS; cell_j++){
            if (cell_i == cell_j) continue;
            Geometry::TetConnectivity tc_j = tet_connect.at(cell_j);

            cell_j_compare = {tc_j.a, tc_j.b, tc_j.c, tc_j.d};
            std::sort(cell_j_compare.begin(), cell_j_compare.end());
            std::set<Index> set_j(cell_j_compare.begin(), cell_j_compare.end());

            vector<Index> common;
            std::set_intersection(set_i.begin(), set_i.end(), set_j.begin(), set_j.end(), std::back_inserter(common));
            if (common.size() >= 3)
            {
                neigbour_found = true;
                break;
            }
        }
        if (!neigbour_found){
            //
            /*for later:
            should have:
            - one function to check if face ij is already created
            - one function to find cell j of face ij given cell i (see above)
             if not it should add ij to the correct boundary patch.
            */
        }

        //
        //Index tet_face_nodes[4][3] = 
        // const vector<vector<Index>> tet_face_nodes = {{tc.a, tc.b, tc.c},  
        //                                             {tc.a, tc.b, tc.d},
        //                                             {tc.a, tc.d, tc.c},
        //                                             {tc.b, tc.c, tc.d}};
       
        // for (const auto& face_ij : tet_face_nodes){
        //     cell_i_compare = face_ij;
        //     std::sort(cell_i_compare.begin(), cell_i_compare.end());

        //     bool neigbour_cell_found = false;
        //     for (Index cell_j=0; cell_j<N_CELLS; cell_j++){
        //         if (cell_j == cell_i) continue;
        //         cell_j_compare = 
        //     }
        // }


    Geometry::Triangle tri_A{a,b,c};
        Geometry::Triangle tri_B{a,b,c};
        Geometry::Triangle tri_C{a,b,c};
        Geometry::Triangle tri_D{a,b,c};
        

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
        cout << "BC_type: "<<BC_type<< "N_S_E"<<N_surface_elements <<endl;
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





    // cout <<"N_N "<<N_NODES<<",N_E "<<N_ELEMENTS<<endl;
    // cout << "nodes:\n";
    // for (auto n : nodes) cout << n <<endl<<endl;;
    // cout <<"tetrahedra:\n";
    // for (auto e : tetrahedra) cout << e.a << ","<<e.b<<","<<e.c<<","<<e.d<<endl;

    // cout << "patches\n";
    // for (auto p : patches){
    //     cout << "bc = "<<(int)p.boundary_type<<endl;
    //     for (auto t : p.triangles){
    //         cout << t.a << ","<<t.b <<","<<t.c<<endl;
    //     }
    // }
    

}