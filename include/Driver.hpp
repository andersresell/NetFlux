#pragma once

#include "Grid.hpp"
#include "Output.hpp"
#include "Solver.hpp"


class Driver{
    
    unique_pointer<geom::Grid> grid;
    unique_pointer<EulerSolver> solver;
    unique_pointer<EulerOutput> output;

public:
    Driver(Config& config);

    void solve(Config& config);


};

