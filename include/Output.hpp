#pragma once
#include "Grid.hpp"
#include "Solver.hpp"

/*This class has been written more generically, but it really wasnt necessary at this point,
only flow::EulerVar is assumed*/

class BaseOutput{
public:
    void write_vtk_ascii(const Config& config, bool write_grid_only=false);
protected:
    virtual void write_vtk_ascii_grid(const Config& config, string filename) = 0;
    virtual void write_vtk_ascii_cell_data(const Config& config, string filename);
    
};

template<ShortIndex N_EQS>
class Output : public BaseOutput{
protected:
    const geom::Grid& grid;
    const FlowField<N_EQS>& solution;

    void write_vtk_ascii_grid(const Config& config, string filename) override;

    Output(const geom::Grid& grid, const EulerField& solution);

public:

};


class EulerOutput : public Output<N_EQS_EULER>{

    void write_vtk_ascii_cell_data(const Config& config, string filename) override;

public:
    EulerOutput(const geom::Grid& grid, const EulerField& solution) : Output(grid, solution) {}


};