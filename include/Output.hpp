#pragma once
#include "geometry/FV_Grid.hpp"
#include "Solver.hpp"

/*This class has been written more generically, but it really wasnt necessary at this point,
only flow::EulerVar is assumed*/

class Output
{
private:
    const geometry::PrimalGrid &primal_grid_glob;
    const geometry::PrimalGrid &primal_grid;
    const Vector<unique_ptr<Solver>> &solvers;
    Vector<unique_ptr<VecField>> consvars_glob;
    void write_vtk_ascii_grid(const Config &config, const string &filename);

public:
    Output(const geometry::PrimalGrid &primal_grid_glob,
           const geometry::PrimalGrid &primal_grid,
           const Vector<unique_ptr<Solver>> &solvers,
           const Config &config);
    void write_vtk_ascii(const Config &config);
    void write_vtk_ascii_debug(const Config &config, const string &filename);
    void gather_cell_data(ShortIndex i_solver);
};

struct EulerOutput
{
    static void write_vtk_ascii_cell_data(const Config &config, const string &filename, const VecField &consvars);
};