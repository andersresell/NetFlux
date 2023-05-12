#include "../include/SolverData.hpp"

EulerSolverData::EulerSolverData(const Config& config){

    Index N_CELLS = config.get_N_TOTAL_CELLS();

    //should specify intial values for solution here
    solution = make_unique<VecField>(N_CELLS, N_EQS_EULER);
    
    //solution_old = make_unique<EulerVecField>();
}
