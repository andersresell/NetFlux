#include "../include/Output.hpp"

Output::Output(const Geom::Grid& grid)
    : grid{grid}
{

}


void Output::write_vtk_ascii(const Config& config){
    //Maybe use data from config to create the filename
    string basename = "output"; //tmp
    string filename = basename + ".vtk"; 
    std::ofstream ost{filename};
    if (!ost){
        std::cerr << "Couldn't open file " << filename << endl;
        exit(1);
    }

    const Index N_NODES = config.get_N_NODES();
    const Index N_TETS = config.get_N_TETS();
    const Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();
    const ShortIndex VTK_TET_TYPE = 10, VTK_TRI_TYPE = 5;
    const auto& nodes = grid.get_nodes();
    const auto& tet_connectivity = grid.get_tet_connect();

    ost << "# vtk DataFile Version 3.0\n"
        << basename << "\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << N_NODES << " double\n";
    
    /*Writing the grid*/
    for (const auto& node : nodes) 
        ost << node.x() << " " << node.y() << " " << node.z() << "\n";

    ost << "\nCELLS " << N_TETS << " " << (1 + N_TET_FACES) * N_TETS << "\n";
    for (const auto& tc : tet_connectivity) 
        ost << N_TET_FACES << " " << tc.a() << " " << tc.b() << " " << tc.c() << " " << tc.d() << "\n";
    ost << "CELL_TYPES " << N_TETS << "\n";
    for (Index i{0}; i<N_TETS; i++)
        ost << VTK_TET_TYPE << "\n";
    
    /*Writing the cell data*/
    //ost << "\n CELL DATA " << N_INTERIOR_CELLS << "\n"; 

}