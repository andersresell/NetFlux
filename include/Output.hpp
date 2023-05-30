#pragma once
#include "Grid.hpp"
#include "Solver.hpp"

/*This class has been written more generically, but it really wasnt necessary at this point,
only flow::EulerVar is assumed*/

class Output{
protected:
    const geom::Grid& grid;
    const Vector<unique_ptr<Solver>>& solvers; 

    virtual void write_vtk_ascii_grid(const Config& config, string filename);
    
public:
    Output(const geom::Grid& grid, const Vector<unique_ptr<Solver>>& solvers);
    void write_vtk_ascii(const Config& config, bool write_grid_only=false);



};



struct EulerOutput {
    static void write_vtk_ascii_cell_data(const Config& config, const string& filename, const VecField& solution) ;
};