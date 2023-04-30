#pragma once
#include "includes.hpp"
#include "Utilities.hpp"

class Config {

    Index N_NODES,
          N_TETS;

    Index N_INTERIOR_CELLS,
          N_TOTAL_CELLS,
          N_FACES;
    bool grid_data_set{false};

    string mesh_filename; 

    map<string, BoundaryType> map_patch_BC; //map from each patch to the bc type applied


public:
    Config(string config_filename); 

    Index get_N_NODES() const {return N_NODES;}
    Index get_N_TETS() const {return N_TETS;}
    Index get_N_INTERIOR_CELLS() const {return N_INTERIOR_CELLS;}
    Index get_N_TOTAL_CELLS() const { return N_TOTAL_CELLS;}
    Index get_N_GHOST_CELLS() const { return get_N_TOTAL_CELLS() - get_N_INTERIOR_CELLS();}
    Index get_N_FACES() const {return N_FACES;}

    void set_grid_data(Index N_NODES, 
                       Index N_INTERIOR_CELLS, 
                       Index N_TOTAL_CELLS, 
                       Index N_FACES);
    
    string get_mesh_filename() const {return mesh_filename;}


    BoundaryType get_boundary_type(string patch_name) const {return map_patch_BC.at(patch_name);}
};

