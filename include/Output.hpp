#pragma once
#include "Grid.hpp"
#include "Solver.hpp"

/*This class has been written more generically, but it really wasnt necessary at this point,
only flow::EulerVar is assumed*/

template<typename FlowVar>
class Output{
protected:
    const geom::Grid& grid;
    const Vector<FlowVar>& solution;

    void write_vtk_ascii_grid(const Config& config, string filename);
    virtual void write_vtk_ascii_cell_data(const Config& config, string filename) = 0;

    Output(const geom::Grid& grid, const Vector<FlowVar>& solution);

public:
    void write_vtk_ascii(const Config& config, bool write_grid_only=false);

};


class EulerOutput : public Output<flow::EulerVar>{

    void write_vtk_ascii_cell_data(const Config& config, string filename) override;

public:
    EulerOutput(const geom::Grid& grid, const Vector<flow::EulerVar>& solution) : Output(grid, solution) {}


};