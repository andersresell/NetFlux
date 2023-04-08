
#include "../include/Utilities.hpp"
#include <iostream>
#include "../include/Grid.hpp"

int main(){
    Config c;
    Grid g{c};
    g.read_mesh("test_mesh.c3d");

    cout << endl << endl;
}