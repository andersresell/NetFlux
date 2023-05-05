#pragma once

#include "Grid.hpp"
#include "Output.hpp"
#include "Solver.hpp"


class Driver{
    
    unique_ptr<geom::Grid> grid;
    Vector<unique_ptr<Solver>> solver;

    unique_ptr<Output> output;

public:
    Driver(Config& config);

    void solve(Config& config);


};

