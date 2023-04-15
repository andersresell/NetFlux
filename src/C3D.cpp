
#include "../include/Utilities.hpp"
#include <iostream>
#include "../include/Grid.hpp"

int main(){

    
    cout << "hello jævel\n\n";
    Config c;

    Grid g{c};  

    g.create_grid(c);

    g.print_grid();
}