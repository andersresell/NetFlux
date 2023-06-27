#include "../include/SolverData.hpp"

SolverData::SolverData(const Config &config, ShortIndex n_eqs)
{
    Index N_TOTAL_CELLS = config.get_N_TOTAL_CELLS();
    Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();

    solution = make_unique<VecField>(N_TOTAL_CELLS, n_eqs);
    solution_old = make_unique<VecField>(N_TOTAL_CELLS, n_eqs);
    primvars = make_unique<VecField>(N_TOTAL_CELLS, n_eqs);
    flux_balance = make_unique<VecField>(N_INTERIOR_CELLS, n_eqs);

    primvars_gradient = make_unique<GradField>(N_INTERIOR_CELLS, n_eqs);
    primvars_limiter = make_unique<VecField>(N_INTERIOR_CELLS, n_eqs);
    primvars_max = make_unique<VecField>(N_INTERIOR_CELLS, n_eqs);
    primvars_min = make_unique<VecField>(N_INTERIOR_CELLS, n_eqs);
}

EulerSolverData::EulerSolverData(const Config &config) : SolverData(config, N_EQS_EULER)
{
    Delta_S.resize(config.get_N_TOTAL_CELLS());

    switch (config.get_initial_cond_option())
    {
    case InitialConditionOption::Freestream:
        set_freestream_values(config);
        break;
    default:
        FAIL_MSG("Error, illegal initial condition option\n");
    }
}

void EulerSolverData::set_primvars(const VecField &cons_vars, const Config &config)
{
    assert(cons_vars.size() == primvars->size() && cons_vars.get_N_EQS() == primvars->get_N_EQS());

    for (Index i{0}; i < config.get_N_INTERIOR_CELLS(); i++)
    {
        const EulerVecMap U = cons_vars.get_variable<EulerVec>(i);
        EulerVecMap V = primvars->get_variable<EulerVec>(i);
        EulerEqs::cons_to_prim(U, V);
    }
}

void EulerSolverData::set_freestream_values(const Config &config)
{
    EulerVec V_inf{config.get_primvars_inf()};
    EulerVec U_inf{};
    Index first{0}, last{config.get_N_INTERIOR_CELLS()};
    EulerEqs::prim_to_cons(V_inf, U_inf);
    primvars->set_constant_field_segment(V_inf, first, last);
    solution->set_constant_field_segment(U_inf, first, last);
}

void ValidityChecker::check_flux_balance_validity(const Config &config, const VecField &flux_balance)
{
    if (!config.check_physical_validity())
        return;

    Index invalid_cells = check_field_validity(flux_balance);

    if (invalid_cells > 0)
    {
        throw std::runtime_error(std::to_string(invalid_cells) + " cells with unphysical \  
        values detected in the flux_balance field of the " +
                                 get_solver_name() + " solver");
    }
}

Index ValidityChecker::check_field_validity(const VecField &field)
{
    Index invalid_cells{0};
    for (Index i{0}; i < field.size(); i++)
    {
        for (ShortIndex j{0}; j < field.get_N_EQS(); j++)
        {
            if (!num_is_valid(field(i, j)))
            {
                invalid_cells += 1;
                break;
            }
        }
    }
    return invalid_cells;
}

Index EulerValidityChecker::check_primvars(const VecField &V)
{
    assert(V.get_N_EQS() == N_EQS_EULER);
    Index invalid_cells{0};
    for (Index i{0}; i < V.size(); i++)
    {
        /*Checking density and pressure*/
        if (!num_is_valid_and_pos(V(i, 0)) || !num_is_valid_and_pos(V(i, 4)))
        {
            invalid_cells += 1;
            continue;
        }
        /*Checking velocity*/
        for (ShortIndex i_dim{1}; i_dim <= N_DIM; i_dim++)
        {
            if (!num_is_valid(V(i, i_dim)))
            {
                invalid_cells += 1;
                break;
            }
        }
    }
    return invalid_cells;
}