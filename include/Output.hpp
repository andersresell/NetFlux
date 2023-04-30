#pragma once
#include "includes.hpp"
#include "Utilities.hpp"
#include "Grid.hpp"

class Output{

    const Geom::Grid& grid;

public:
    Output(const Geom::Grid& grid);


    void write_vtk_ascii(const Config& config);
    
};
