#include "../include/SolverData.hpp"

EulerSolverData::EulerSolverData(const Config& config){

    Index N_CELLS = config.get_N_TOTAL_CELLS();

    //should specify intial values for solution here
    solution = make_unique<VecField>(N_CELLS, N_EQS_EULER);
    
    Delta_S.resize(N_CELLS);
    //solution_old = make_unique<EulerVecField>();
}



void EulerSolverData::set_primvars(const VecField& cons_vars, const Config& config){
    
    assert(cons_vars.size() == primvars->size() && cons_vars.get_N_EQS() == primvars->get_N_EQS());

    EulerVecMap U{nullptr}, V{nullptr};
    for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
        U = cons_vars.get_variable<EulerVec>(i);
        V = primvars->get_variable<EulerVec>(i);
        EulerEqs::cons_to_prim(U, V);
    }
}