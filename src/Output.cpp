
#include "../include/Output.hpp"





Output::Output(const geom::Grid& grid, const Vector<const Solver&> solvers)
    : grid{grid}, solvers{move(solvers)}
{

}

void Output::write_vtk_ascii(const Config& config, bool write_grid_only){
    const string& filename = config.get_unsteady_vtk_filename();

    write_vtk_ascii_grid(config, filename);
    if (!write_grid_only) {
        for (const Solver& solver : solvers){
            
            switch (solver.get_solver_type()){
                case SolverType::Euler:{
                    auto solution = static_cast<const EulerVecField&>(solver.get_solver_data().get_solution());
                    EulerOutput::write_vtk_ascii_cell_data(config, filename, solution);
                    break;
                }
                default:
                    FAIL_MSG("Error: Illegal SolverType\n");
            }
        }
        
        
    }
}


void Output::write_vtk_ascii_grid(const Config& config, string filename){
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


void EulerOutput::write_vtk_ascii_cell_data(const Config& config, string filename, const EulerVecField& solution){
    
    std::ofstream ost{filename, std::ios::app};
    FAIL_IF_MSG(!ost, "Couldn't open file " + filename);

    const Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();
    const auto& U = solution.cell_values;

    ost << "\n CELL DATA " << N_INTERIOR_CELLS << "\n"; 
    
    ost << "\nSCALARS density double 1\n";
    for (Index i{0}; i< N_INTERIOR_CELLS; i++)
        ost << U[i][0] << "\n";

    ost << "\nVECTORS velocity double\n"; 
    for (Index i{0}; i< N_INTERIOR_CELLS; i++){
        double density = U[i][0];
        ost << U[i][1]/density << " " << U[i][2]/density << " " << U[i][3]/density << "\n";
    }

    ost << "\nSCALARS pressure double 1\n";
    for (Index i{0}; i< N_INTERIOR_CELLS; i++)
        ost << EulerField::pressure(U[i]) << "\n";

}