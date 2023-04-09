#pragma once
#include "includes.hpp"
#include "Utilities.hpp"


class Config {

    Index N_INTERIOR_CELLS,
          N_TOTAL_CELLS,
          N_FACES;
    bool faces_set = false, interior_cells_set = false, total_cells_set;
    string mesh_filename; 
public:
    Index get_N_INTERIOR_CELLS() const {assert(interior_cells_set); return N_INTERIOR_CELLS;}
    Index get_N_TOTAL_CELLS() const {assert(total_cells_set); return N_TOTAL_CELLS;}
    Index get_N_GHOST_CELLS() const {return get_N_TOTAL_CELLS() - get_N_INTERIOR_CELLS();}

    Index get_N_FACES() const {return N_FACES;}
    string get_mesh_filename() const {return mesh_filename;}
    


    void set_N_INTERIOR_CELLS(Index val) {assert(!interior_cells_set); N_INTERIOR_CELLS = val; interior_cells_set = true;}
    void set_N_TOTAL_CELLS(Index val) {assert(!total_cells_set); N_TOTAL_CELLS = val; total_cells_set = true;}
    void set_N_FACES(Index val) {assert(!faces_set); N_FACES = val; faces_set = true;}

};

