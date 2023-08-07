#pragma once

#include "geometry/PrimalGrid.hpp"
#include "geometry/FV_Grid.hpp"
#include "Output.hpp"
#include "Solver.hpp"

class Driver
{

    Config &config;

    unique_ptr<geometry::PrimalGrid> primal_grid_glob;

    unique_ptr<geometry::PrimalGrid> primal_grid;

    unique_ptr<geometry::FV_Grid> FV_grid;

    vector<unique_ptr<Solver>> solvers;

    unique_ptr<Output> output;

public:
    Driver(Config &config);

    void solve();
};
