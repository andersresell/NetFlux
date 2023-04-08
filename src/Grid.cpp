#include "../include/Grid.hpp"

Grid::Grid(Config config){

}

void Grid::read_mesh(string mesh_file){
   
    std::ifstream ist{mesh_file};
    std::stringstream ss;

    if (!ist) {
        std::cerr << "Error: couldn't open the mesh file: "+mesh_file;
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

    vector<Vec3> nodes; nodes.resize(N_NODES);
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
    
    struct Tetrahedron{
        ShortIndex a,b,c,d;
    };
    vector<Tetrahedron> tetrahedra; tetrahedra.resize(N_ELEMENTS);
    while (getline(ist, line)) {
        if (line.size() != 0){
            assert(line == "#tetrahedra");
            break;
        } 
    }
    //reading elements (assuming only tetrahedra for now)
    for (Index i{0}; i<N_ELEMENTS; i++){
        Tetrahedron t;
        ist >> t.a >> t.b >> t.c >> t.d;
        tetrahedra.at(i) = t;
    }

    struct Triangle{
        ShortIndex a,b,c;
    };
    struct Patch{
        BoundaryType boundary_type;
        vector<Triangle> triangles;
    };

    while (getline(ist, line)) {
        if (line.size() != 0){
            assert(line == "#patches");
            break;
        } 
    }
    string BC_type;
    Index N_surface_elements;    
    vector<Patch> patches; //Might consider preallocating memory if this one is slow
    //Reading boundary patches
    while (ist >> BC_type >> N_surface_elements){
        cout << "BC_type: "<<BC_type<< "N_S_E"<<N_surface_elements <<endl;
        BoundaryType boundary_type = boundary_type_from_string.at(BC_type);
        
        Triangle tri;
        Patch p;
        p.boundary_type = boundary_type; 
        p.triangles.resize(N_surface_elements);
        for (Index i{0}; i<N_surface_elements; i++){
            ist >> tri.a >> tri.b >> tri.c;
            p.triangles.at(i) = tri;
        }
        patches.emplace_back(p);
    }





    cout <<"N_N "<<N_NODES<<",N_E "<<N_ELEMENTS<<endl;
    cout << "nodes:\n";
    for (auto n : nodes) cout << n <<endl<<endl;;
    cout <<"tetrahedra:\n";
    for (auto e : tetrahedra) cout << e.a << ","<<e.b<<","<<e.c<<","<<e.d<<endl;

    cout << "patches\n";
    for (auto p : patches){
        cout << "bc = "<<(int)p.boundary_type<<endl;
        for (auto t : p.triangles){
            cout << t.a << ","<<t.b <<","<<t.c<<endl;
        }
    }
    

}