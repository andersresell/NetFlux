#include "../include/SolverData.hpp"


SolverData::SolverData(const Config& config, ShortIndex n_eqs){
    Index N_TOTAL_CELLS = config.get_N_TOTAL_CELLS();
    Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();
    
    solution = make_unique<VecField>(N_TOTAL_CELLS, n_eqs);
    solution_old = make_unique<VecField>(N_TOTAL_CELLS, n_eqs);
    primvars = make_unique<VecField>(N_TOTAL_CELLS, n_eqs);
    flux_balance = make_unique<VecField>(N_INTERIOR_CELLS, n_eqs);
   

    //solution_old = make_unique<EulerVecField>();

}

EulerSolverData::EulerSolverData(const Config& config) : SolverData(config, N_EQS_EULER){
    Delta_S.resize(config.get_N_TOTAL_CELLS());

    switch (config.get_initial_cond_option()){
        case InitialConditionOption::Freestream:
            set_freestream_values(config);
            break;
        default:
            FAIL_MSG("Error, illegal initial condition option\n");
    }
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


void EulerSolverData::set_freestream_values(const Config& config){
    EulerVec V_inf{config.get_primvars_inf()};
    EulerVec U_inf{};
    Index first{0}, last{config.get_N_INTERIOR_CELLS()};
    EulerEqs::prim_to_cons(V_inf, U_inf);
    primvars->set_constant_field_segment(V_inf, first, last);
    solution->set_constant_field_segment(U_inf, first, last);
}