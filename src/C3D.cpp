
#include "../include/Utilities.hpp"
#include <iostream>
#include "../include/Grid.hpp"

int main(){


    cout << "hello jævel\n\n";

    // cout << "N_DIM " << N_DIM<<endl;
    // cout << "N_TET_NODES "<<N_TET_NODES<<endl;
    // cout << "N_TET_FACES "<<N_TET_FACES<<endl;

    Config c{"/home/anders/dev/Compress3D/test_mesh.c3d"};

    Grid g{c};  

    g.create_grid(c);

    //g.print_grid();
}