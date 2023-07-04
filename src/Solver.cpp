#include "../include/Solver.hpp"

using namespace geom;

Solver::Solver(const geom::Grid &grid, const Config &config)
    : grid{grid}
{
    create_BC_container(config);
}

void Solver::create_BC_container(const Config &config)
{
    for (const auto &patch : grid.get_patches())
    {
        unique_ptr<BoundaryCondition> BC;
        switch (patch.boundary_type)
        {
        case BoundaryType::NoSlipWall:
            BC = make_unique<BC_NoSlipWall>();
            break;
        case BoundaryType::SlipWall:
            BC = make_unique<BC_SlipWall>();
            break;
        case BoundaryType::FarField:
            BC = make_unique<BC_FarField>(config);
            break;
        default:
            assert(false);
        }
        BC_container.emplace_back(move(BC));
    }
}

void Solver::step(const Config &config)
{

    assert(config.get_time_integration_type() == TimeIntegrationType::Explicit); // Remove if implementing implicit

    switch (config.get_time_scheme())
    {
    case TimeScheme::ExplicitEuler:
        explicit_euler(config);
        break;
    case TimeScheme::TVD_RK3:
        TVD_RK3(config);
        break;
    default:
        assert(false);
    }

    solver_data->get_solution_old() = solver_data->get_solution();
}

void Solver::evaluate_flux_balance(const Config &config, const VecField &cons_vars)
{
    assert(validity_checker->valid_consvars_interior(cons_vars));

    solver_data->get_flux_balance().set_zero();

    solver_data->set_primvars(cons_vars, config);

    assert(validity_checker->valid_primvars_interior(solver_data->get_primvars()));
    set_constant_ghost_values(config);
    assert(validity_checker->valid_primvars_ghost(solver_data->get_primvars()));

    validity_checker->write_debug_info(solver_data->get_primvars(), "Primvars");

    if (config.get_spatial_order() == SpatialOrder::Second)
    {
        evaluate_gradient(config);

        if (config.get_limiter() != Limiter::NONE)
            evaluate_limiter(config);
    }

    evaluate_inviscid_fluxes(config);

    evaluate_viscous_fluxes(config);

    validity_checker->check_flux_balance_validity(config, solver_data->get_flux_balance());
}

void Solver::explicit_euler(const Config &config)
{
    const Scalar dt = config.get_delta_time();
    VecField &U = solver_data->get_solution();
    VecField &R = solver_data->get_flux_balance();
    const ShortIndex N_EQS = solver_data->get_N_EQS();
    const auto &cells = grid.get_cells();
    Index i, j;

    assert(U.get_N_EQS() == N_EQS && R.get_N_EQS() == N_EQS);
    assert(U.size() == config.get_N_INTERIOR_CELLS() && R.size() == config.get_N_INTERIOR_CELLS());

    /*--------------------------------------------------------------------
     U_n+1 = U_n + dt /Omega * R(U_n)
    --------------------------------------------------------------------*/

    evaluate_flux_balance(config, U);
    for_all(U, i, j)
        U(i, j) += dt / cells[i].cell_volume * R(i, j);
}

void Solver::TVD_RK3(const Config &config)
{
    const Scalar dt = config.get_delta_time();
    VecField &U = solver_data->get_solution();
    VecField &U_old = solver_data->get_solution_old();
    VecField &R = solver_data->get_flux_balance();
    const ShortIndex N_EQS = solver_data->get_N_EQS();
    const auto &cells = grid.get_cells();
    Index i, j;

    assert(U.get_N_EQS() == N_EQS && U_old.get_N_EQS() == N_EQS && R.get_N_EQS() == N_EQS);
    assert(U.size() == config.get_N_INTERIOR_CELLS() && U_old.size() == config.get_N_INTERIOR_CELLS());
    assert(R.size() == config.get_N_INTERIOR_CELLS());

    /*--------------------------------------------------------------------
    U_1 = U_n + dt /Omega * R(U_n)
    U_2 = 3/4 * U_n + 1/4 *U_1 + 1/4 * dt / Omega * R(U_1)
    U_n+1 = 1/3 * U_n + 2/3 * U_2 + 2/3 * dt / Omega * R(U_2)
    --------------------------------------------------------------------*/

    evaluate_flux_balance(config, U);
    for_all(U, i, j)
        U(i, j) += dt / cells[i].cell_volume * R(i, j);

    evaluate_flux_balance(config, U);
    for_all(U, i, j)
        U(i, j) = 3.0 / 4.0 * U_old(i, j) + 1.0 / 4.0 * U(i, j) + 1.0 / 4.0 * dt / cells[i].cell_volume * R(i, j);

    evaluate_flux_balance(config, U);
    for_all(U, i, j)
        U(i, j) = 1.0 / 3.0 * U_old(i, j) + 2.0 / 3.0 * U(i, j) + 2.0 / 3.0 * dt / cells[i].cell_volume * R(i, j);
}

EulerSolver::EulerSolver(const Config &config, const geom::Grid &grid) : Solver(grid, config)
{
    solver_data = make_unique<EulerSolverData>(config);
    validity_checker = make_unique<EulerValidityChecker>(config);

    calc_Delta_S(config);
}

void EulerSolver::evaluate_inviscid_fluxes(const Config &config)
{
    VecField &flux_balance = solver_data->get_flux_balance();

    // const VecField& primvars = solver_data->get_primvars();
    // const GradField& primvars_grad = solver_data->get_primvars_gradient();
    // const VecField& primvars_limiter = solver_data->get_primvars_limiter();

    Index N_INTERIOR_FACES = config.get_N_INTERIOR_FACES();
    SpatialOrder spatial_order = config.get_spatial_order();

    Index i, j;
    const auto &faces = grid.get_faces();
    const auto &patches = grid.get_patches();

    // InvFluxFunction inv_flux_func = NumericalFlux::get_inviscid_flux_function(config);

    EulerSolverData &euler_data = dynamic_cast<EulerSolverData &>(*solver_data);
    EulerVecMap U_L = euler_data.get_U_L_map();
    EulerVecMap U_R = euler_data.get_U_R_map();
    EulerVecMap V_L = euler_data.get_V_L_map();
    EulerVecMap V_R = euler_data.get_V_R_map();
    EulerVecMap Flux_inv = euler_data.get_Flux_inv_map();

    NumericalFlux::InvFluxFunction numerical_flux_func = NumericalFlux::get_inviscid_flux_function(config.get_inv_flux_scheme());

    /*First interior cells*/
    for (Index ij{0}; ij < N_INTERIOR_FACES; ij++)
    {
        i = faces[ij].i;
        j = faces[ij].j;
        const Vec3 &S_ij = faces[ij].S_ij;
        const Vec3 &r_im = faces[ij].r_im;
        const Vec3 &r_jm = faces[ij].r_jm;

        calc_reconstructed_value(i, V_L, r_im, spatial_order);
        calc_reconstructed_value(j, V_R, r_jm, spatial_order);

        EulerEqs::prim_to_cons(V_L, U_L);
        EulerEqs::prim_to_cons(V_R, U_R);

        numerical_flux_func(U_L, U_R, S_ij, Flux_inv);
        // cout << "F_inf " << endl
        //      << Flux_inv << endl
        //      << endl;

        flux_balance.get_variable<EulerVec>(i) -= Flux_inv;
        flux_balance.get_variable<EulerVec>(j) += Flux_inv;
    }

    /*Then boundaries. Here ghost cells has to be assigned based on the boundary conditions.
    This is handled patch-wise*/

    for (Index i_patch{0}; i_patch < patches.size(); i_patch++)
    {
        const auto &patch = patches[i_patch];
        auto &boundary_condition = BC_container[i_patch];
        // BoundaryCondition::BC_function BC_func = BoundaryCondition::get_BC_function(patch.boundary_type);

        for (Index ij{patch.FIRST_FACE}; ij < patch.FIRST_FACE + patch.N_FACES; ij++)
        {

            i = faces[ij].i;
            j = faces[ij].j;
            const Vec3 &S_ij = faces[ij].S_ij;
            const Vec3 &r_im = faces[ij].r_im;

            calc_reconstructed_value(i, V_L, r_im, spatial_order);

            boundary_condition->calc_ghost_val(V_L, V_R, S_ij);

            EulerEqs::prim_to_cons(V_L, U_L);
            EulerEqs::prim_to_cons(V_R, U_R);

            numerical_flux_func(U_L, U_R, S_ij, Flux_inv);

            assert(validity_checker->valid_boundary_flux(Flux_inv.data(), patch.boundary_type));

            flux_balance.get_variable<EulerVec>(i) -= Flux_inv;
        }
    }
}

void EulerSolver::calc_timestep(Config &config)
{
    // --------------------------------------------------------------------
    // Implementing Method 2 in "Time Step on Unstructured Grids" in Blazek
    // --------------------------------------------------------------------

    /*Delta S only needs updating when the grid is moved*/
    if (config.get_grid_motion())
        calc_Delta_S(config);

    const EulerSolverData &euler_data = dynamic_cast<const EulerSolverData &>(*solver_data);
    const Vector<Vec3> &Delta_S = euler_data.get_Delta_S();

    const Scalar CFL = config.get_CFL();
    auto cells = grid.get_cells();

    const auto &primvars = solver_data->get_primvars();

    Scalar c, volume, spec_rad_x, spec_rad_y, spec_rad_z;

    Scalar delta_time = std::numeric_limits<Scalar>::max(); // Large number

    for (Index i{0}; i < config.get_N_INTERIOR_CELLS(); i++)
    {
        const EulerVecMap V = primvars.get_variable<EulerVec>(i);
        c = EulerEqs::sound_speed_primitive(V);

        volume = cells[i].cell_volume;

        spec_rad_x = (abs(V[1]) + c) * Delta_S[i].x();
        spec_rad_y = (abs(V[2]) + c) * Delta_S[i].y();
        spec_rad_z = (abs(V[3]) + c) * Delta_S[i].z();

        delta_time = std::min(delta_time, CFL * volume / (spec_rad_x + spec_rad_y + spec_rad_z));
    }
    if (!num_is_valid_and_pos(delta_time))
        throw std::runtime_error("Invalid dt calculated (dt = " + std::to_string(delta_time) + ")");

    if (delta_time > 0.01)
        cout << "Warning: Delta time is high (dt = " << delta_time << ")\n";

    config.set_delta_time(delta_time);
}

void EulerSolver::calc_Delta_S(const Config &config)
{

    const auto &faces = grid.get_faces();
    EulerSolverData &euler_data = dynamic_cast<EulerSolverData &>(*solver_data);

    Vector<Vec3> &Delta_S = euler_data.get_Delta_S();
    std::fill(Delta_S.begin(), Delta_S.end(), Vec3::Zero());
    Index i, j;
    Vec3 tmp;

    for (Index ij{0}; ij < config.get_N_INTERIOR_FACES(); ij++)
    {
        i = faces[ij].i;
        j = faces[ij].j;

        tmp = 0.5 * faces[ij].S_ij.cwiseAbs();

        Delta_S[i] += tmp;
        Delta_S[j] += tmp;
    }
}

void EulerSolver::set_constant_ghost_values(const Config &config)
{
    const auto &patches = grid.get_patches();
    const auto &faces = grid.get_faces();
    VecField &primvars = solver_data->get_primvars();
    Index i, j;

    // for (const auto &patch : patches)
    for (Index i_patch{0}; i_patch < patches.size(); i_patch++)
    {
        const auto &patch = patches[i_patch];
        auto &boundary_condition = BC_container[i_patch];

        for (Index ij{patch.FIRST_FACE}; ij < patch.FIRST_FACE + patch.N_FACES; ij++)
        {
            i = faces[ij].i;
            j = faces[ij].j;
            assert(i < config.get_N_INTERIOR_CELLS() && j >= config.get_N_INTERIOR_CELLS());

            const Vec3 &S_ij = faces[ij].S_ij;
            // primvars.map_to_variable<EulerVec>(j) =
            //     BoundaryCondition::calc_ghost_val<BC_type, EulerVec>(primvars.map_to_variable<EulerVec>(i), S_ij);

            // Test if this works
            EulerVecMap V_i = primvars.get_variable<EulerVec>(i);
            EulerVecMap V_j = primvars.get_variable<EulerVec>(j);

            boundary_condition->calc_ghost_val(V_i, V_j, S_ij);
        }
    }
}

void EulerSolver::evaluate_gradient(const Config &config)
{
    const VecField &primvars = solver_data->get_primvars();
    GradField &primvars_grad = solver_data->get_primvars_gradient();

    switch (config.get_grad_scheme())
    {
    case GradientScheme::GreenGauss:
        Gradient::calc_green_gauss_gradient<N_EQS_EULER>(config, grid, primvars, primvars_grad);
        break;
    default:
        FAIL_MSG("Selected gradient scheme not implemented\n");
    }
}

void EulerSolver::evaluate_limiter(const Config &config)
{

    const VecField &primvars = solver_data->get_primvars();
    const GradField &primvars_grad = solver_data->get_primvars_gradient();
    VecField &primvars_limiter = solver_data->get_primvars_limiter();
    VecField &primvars_max = solver_data->get_primvars_max();
    VecField &primvars_min = solver_data->get_primvars_min();

    switch (config.get_limiter())
    {
    case Limiter::Barth:
        Reconstruction::calc_max_and_min_values<N_EQS_EULER>(config,
                                                             grid,
                                                             primvars,
                                                             primvars_max,
                                                             primvars_min);
        Reconstruction::calc_barth_limiter<N_EQS_EULER>(config,
                                                        grid,
                                                        primvars,
                                                        primvars_grad,
                                                        primvars_max,
                                                        primvars_min,
                                                        primvars_limiter);
        break;
    default:
        FAIL_MSG("Selected limiter not implemented\n");
    }
}
