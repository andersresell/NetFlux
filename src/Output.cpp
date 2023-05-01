#include "../include/Output.hpp"

template<typename FlowVar>
Output<FlowVar>::Output(const geom::Grid& grid, const Vector<FlowVar>& solution)
    : grid{grid}, solution{solution}
{

}

template<typename FlowVar>
void Output<FlowVar>::write_vtk_ascii_grid(const Config& config, string filename){
    std::ofstream ost{filename};
    FAIL_IF_MSG(!ost, "Couldn't open file " + filename);

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
}


template<typename FlowVar>
void Output<FlowVar>::write_vtk_ascii(const Config& config, bool write_grid_only){
    string filename = config.get_unsteady_vtk_filename();

    write_vtk_ascii_grid(config, filename);
    if (!write_grid_only) write_vtk_ascii_cell_data(config, filename);    
}





void EulerOutput::write_vtk_ascii_cell_data(const Config& config, string filename){

    std::ofstream ost{filename, std::ios::app};
    FAIL_IF_MSG(!ost, "Couldn't open file " + filename);

    const Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();

    ost << "\n CELL DATA " << N_INTERIOR_CELLS << "\n"; 
    
    ost << "\nSCALARS density double 1\n";
    for (Index i{0}; i< N_INTERIOR_CELLS; i++)
        ost << solution[i][0] << "\n";

    ost << "\nVECTORS velocity double\n"; 
    for (Index i{0}; i< N_INTERIOR_CELLS; i++){
        double density = solution[i][0];
        ost << solution[i][1]/density << " " << solution[i][2]/density << " " << solution[i][3]/density << "\n";
    }

    ost << "\nSCALARS pressure double 1\n";
    for (Index i{0}; i< N_INTERIOR_CELLS; i++)
        ost << flow::EulerVar::pressure(solution[i]) << "\n";

}