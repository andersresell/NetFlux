
#include "../include/Output.hpp"





Output::Output(const geom::Grid& grid, const Vector<unique_ptr<Solver>>& solvers)
    : grid{grid}, solvers{solvers}
{

}

void Output::write_vtk_ascii(const Config& config, bool write_grid_only){
    const string& filename = config.get_unsteady_vtk_filename();

    write_vtk_ascii_grid(config, filename);
    if (!write_grid_only) {

        for (const auto& solver : solvers){
            
            const VecField& primvars = solver->get_solver_data().get_primvars();

            switch (solver->get_solver_type()){
                case SolverType::Euler:{
                    EulerOutput::write_vtk_ascii_cell_data(config, filename, primvars);
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
    const ShortIndex VTK_TET_TYPE = 10; //VTK_TRI_TYPE = 5;
    const auto& nodes = grid.get_nodes();
    const auto& tet_connectivity = grid.get_tet_connect();

    ost << "# vtk DataFile Version 3.0\n"
        << "Compress 3D\n"
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


void EulerOutput::write_vtk_ascii_cell_data(const Config& config, const string& filename, const VecField& primvars){
    
    assert(primvars.get_N_EQS() == N_EQS_EULER);

    std::ofstream ost{filename, std::ios::app};
    FAIL_IF_MSG(!ost, "Couldn't open file " + filename);

    const Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();

    ost << "\n CELL DATA " << N_INTERIOR_CELLS << "\n"; 
    
    ost << "\nSCALARS density double 1\n";
    for (Index i{0}; i< N_INTERIOR_CELLS; i++)
        ost << primvars(i,0) << "\n";

    ost << "\nVECTORS velocity double\n"; 
    for (Index i{0}; i< N_INTERIOR_CELLS; i++){
        ost << primvars(i,1) << " " << primvars(i,2) << " " << primvars(i,3) << "\n";
    }

    ost << "\nSCALARS pressure double 1\n";
    for (Index i{0}; i< N_INTERIOR_CELLS; i++)
        ost << primvars(i,4) << "\n";

}