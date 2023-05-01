#include "../include/Solver.hpp"

using namespace flow;

EulerSolver::EulerSolver(const Config& config, const geom::Grid& grid)
    : grid{grid}
{
    //Specify intial val rather
    solution.resize(config.get_N_INTERIOR_CELLS());
    solution_old.resize(config.get_N_INTERIOR_CELLS());
    
}




double EulerSolver::calc_timestep(Config& config){
    double CFL = config.get_CFL();
    for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
        //solution[i].sound_speed();    
    }
}
