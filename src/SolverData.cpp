#include "../include/SolverData.hpp"

SolverData::SolverData(const Config &config, ShortIndex n_eqs)
{
    Index N_CELLS_TOT = config.get_N_CELLS_TOT();
    Index N_CELLS_DOMAIN = config.get_N_CELLS_DOMAIN();
    Index N_CELLS_INT = config.get_N_CELLS_INT();

    solution = make_unique<VecField>(N_CELLS_INT, n_eqs);
    solution_old = make_unique<VecField>(N_CELLS_INT, n_eqs);
    primvars = make_unique<VecField>(N_CELLS_TOT, n_eqs);
    flux_balance = make_unique<VecField>(N_CELLS_INT, n_eqs);

    if (config.get_spatial_order() == SpatialOrder::Second)
    {
        primvars_gradient = make_unique<GradField>(N_CELLS_DOMAIN, n_eqs);
        primvars_limiter = make_unique<VecField>(N_CELLS_DOMAIN, n_eqs);
        *primvars_limiter = 1.0;
        primvars_max = make_unique<VecField>(N_CELLS_DOMAIN, n_eqs); //?correct size
        primvars_min = make_unique<VecField>(N_CELLS_DOMAIN, n_eqs); //?correct size

        n_vecfields_sendrecv_max = 4;  /*primvars, primvars limiter, primvars max+min*/
        n_gradfields_sendrecv_max = 1; /*primvars grad*/
    }
    else
    {
        n_vecfields_sendrecv_max = 1; /*primvars*/
        n_gradfields_sendrecv_max = 0;
    }
}

EulerSolverData::EulerSolverData(const Config &config) : SolverData(config, N_EQS_EULER)
{
    Delta_S.resize(config.get_N_CELLS_TOT());

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
    assert(cons_vars.size() == config.get_N_INTERIOR_CELLS() && primvars->size() == config.get_N_TOTAL_CELLS());
    assert(cons_vars.get_N_EQS() == primvars->get_N_EQS());

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

ValidityChecker::ValidityChecker(const Config &config) : config{config}
{
    std::ofstream ost{DEBUG_LOG_FILE};
    FAIL_IF_MSG(!ost, "Couldn't open file debugging logging file " + string(DEBUG_LOG_FILE));
}

void ValidityChecker::check_flux_balance_validity(const Config &config, const VecField &flux_balance) const
{
    if (!config.check_physical_validity())
        return;

    Index invalid_cells = check_field_validity(flux_balance, 0, config.get_N_INTERIOR_CELLS());

    if (invalid_cells > 0)
    {
        throw std::runtime_error(std::to_string(invalid_cells) + " cells with unphysical " +
                                 "values detected in the flux_balance field of the " +
                                 get_solver_name() + " solver");
    }
}

Index ValidityChecker::check_field_validity(const VecField &field, Index first, Index last) const
{

    assert(first >= 0 && last <= field.size());
    Index invalid_cells{0};
    for (Index i{first}; i < last; i++)
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

Index EulerValidityChecker::check_primvars(const VecField &V, Index first, Index last) const
{

    assert(first >= 0 && last <= V.size());
    assert(V.get_N_EQS() == N_EQS_EULER);
    Index invalid_cells{0};
    for (Index i{first}; i < last; i++)
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

Index EulerValidityChecker::check_consvars(const VecField &U, Index first, Index last) const
{
    assert(first >= 0 && last <= U.size());
    assert(U.get_N_EQS() == N_EQS_EULER);
    Index invalid_cells{0};
    for (Index i{first}; i < last; i++)
    {
        /*Checking density and total energy*/
        if (!num_is_valid_and_pos(U(i, 0)) || !num_is_valid_and_pos(U(i, 4)))
        {
            invalid_cells += 1;
            continue;
        }
        /*Checking momentum*/
        for (ShortIndex i_dim{1}; i_dim <= N_DIM; i_dim++)
        {
            if (!num_is_valid(U(i, i_dim)))
            {
                invalid_cells += 1;
                break;
            }
        }
    }
    return invalid_cells;
}

bool ValidityChecker::valid_primvars_interior(const VecField &V) const
{
    if (check_primvars(V, 0, config.get_N_INTERIOR_CELLS()) > 0)
    {
        write_debug_info(V);
        return false;
    };
    return true;
}

bool ValidityChecker::valid_consvars_interior(const VecField &U) const
{
    if (check_consvars(U, 0, config.get_N_INTERIOR_CELLS()) > 0)
    {
        write_debug_info(U);
        return false;
    };
    return true;
}

bool ValidityChecker::valid_primvars_ghost(const VecField &V) const
{
    if (check_consvars(V, config.get_N_INTERIOR_CELLS(), V.size()) > 0)
    {
        write_debug_info(V);
        return false;
    };
    return true;
}

bool ValidityChecker::valid_consvars_ghost(const VecField &U) const
{
    if (check_consvars(U, config.get_N_INTERIOR_CELLS(), U.size()) > 0)
    {
        write_debug_info(U);
        return false;
    };
    return true;
}

bool ValidityChecker::valid_flux_balance(const VecField &R) const
{
    if (check_field_validity(R, 0, config.get_N_INTERIOR_CELLS()) > 0)
    {
        write_debug_info(R);
        return false;
    };
    return true;
}

void ValidityChecker::write_debug_info(const VecField &U, string name) const
{
    std::ofstream ost{DEBUG_LOG_FILE, std::ios::app};
    FAIL_IF_MSG(!ost, "Couldn't open file debugging logging file " + string(DEBUG_LOG_FILE));

    ost << "\n\nDisplaying vector field " + name + "\n\n";
    for (Index i{0}; i < U.size(); i++)
    {
        for (ShortIndex j{0}; j < U.get_N_EQS(); j++)
        {
            ost << U(i, j);
            if (j < U.get_N_EQS() - 1)
                ost << ", ";
        }
        if (i >= config.get_N_INTERIOR_CELLS())
            ost << " -G";
        ost << "\n";
    }
}

bool EulerValidityChecker::valid_boundary_flux(const Scalar *flux_vals, BoundaryType bc_type) const
{
    const EulerVec flux = EulerVec{flux_vals};
    assert(flux.allFinite() && !flux.hasNaN());

    if (bc_type == BoundaryType::NoSlipWall || bc_type == BoundaryType::SlipWall)
    {
        /*Check that the density and energy flux is zero due to zero normal velocity*/
        if (!is_approx_zero(flux[0]) && !is_approx_zero(flux[4]))
            return false;
    }
    return true;
}