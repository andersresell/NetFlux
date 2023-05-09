#include "../include/SolverData.hpp"

EulerSolverData::EulerSolverData(const Config& config){


    //should specify intial values for solution here
    solution = make_unique<EulerVecField>();
    
    solution_old = make_unique<EulerVecField>();
}
