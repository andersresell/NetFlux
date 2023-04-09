#pragma once

#include "includes.hpp"
#include "Utilities.hpp"
#include "Config.hpp"



class Grid{
    CellContainer cells;
    vector<Patch> patches;
    FaceContainer faces;




public:    
    
    Grid(const Config& config);
    void create_grid(const Config& config);

    void read_mesh(string mesh_filename, 
                   vector<Vec3>& nodes, 
                   vector<Geometry::TetConnectivity>& tet_connect, 
                   vector<Geometry::TriPatchConnectivity>& tri_patch_connect) const;


};