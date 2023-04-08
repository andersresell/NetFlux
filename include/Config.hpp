#pragma once
#include "includes.hpp"
#include "Utilities.hpp"


class Config {

    Index N_CELLS,
          N_FACES;
    bool faces_set = false, cells_set = false;

public:

    void set_N_CELLS(Index N_CELLS) {assert(!cells_set); N_CELLS = N_CELLS; cells_set = true;}
    void set_N_FACES(Index N_FACES) {assert(!faces_set); N_FACES = N_FACES; faces_set = true;}
};

