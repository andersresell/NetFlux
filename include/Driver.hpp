#pragma once

#include "Grid.hpp"
#include "Output.hpp"
#include "Solver.hpp"


class Driver{
    
    unique_pointer<geom::Grid> grid;
    unique_pointer<BaseSolver> solver;
    unique_pointer<BaseOutput> output;

public:
    Driver(Config& config);

    void solve(Config& config);


};

