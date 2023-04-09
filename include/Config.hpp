#pragma once
#include "includes.hpp"
#include "Utilities.hpp"


class Config {

    Index N_CELLS,
          N_FACES;
    bool faces_set = false, cells_set = false;
    string mesh_filename; 
public:
    Index get_N_CELLS() const {return N_CELLS;}
    Index get_N_FACES() const {return N_FACES;}
    string get_mesh_filename() const {return mesh_filename;}
    


    void set_N_CELLS(Index N_CELLS) {assert(!cells_set); N_CELLS = N_CELLS; cells_set = true;}
    void set_N_FACES(Index N_FACES) {assert(!faces_set); N_FACES = N_FACES; faces_set = true;}
};

