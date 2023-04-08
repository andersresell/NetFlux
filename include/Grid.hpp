#pragma once

#include "includes.hpp"
#include "Utilities.hpp"
#include "Config.hpp"



class Grid{
    CellContainer internal_cells,
                  ghost_cells;
    FaceContainer faces;

    Grid(Config config);

    void read_mesh(string mesh_file);
};