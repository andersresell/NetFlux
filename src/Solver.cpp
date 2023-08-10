#include "../include/Solver.hpp"

using namespace geometry;

Solver::Solver(const Config &config, const PrimalGrid &primal_grid, const FV_Grid &FV_grid, PartitionComm &part_comm)
    : primal_grid{primal_grid}, FV_grid{FV_grid}, part_comm{part_comm}
{
    create_BC_container(config);
}

void Solver::create_BC_container(const Config &config)
{
    for (const auto &p : FV_grid.get_patches_boundary())
    {
        unique_ptr<BoundaryCondition> BC;
        switch (p.boundary_type)
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

    if (config.get_spatial_order() == SpatialOrder::First)
    {
        communicate_primvars();
    }
    else
    {
        assert(config.get_spatial_order() == SpatialOrder::Second);

        communicate_primvars_and_gradient();

        evaluate_gradient(config);

        if (config.get_limiter() != Limiter::NONE)
        {
            evaluate_limiter(config);
        }
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
    const auto &cells = FV_grid.get_cells();
    Index i, n_eq;

    assert(U.get_N_EQS() == solver_data->get_N_EQS() && R.get_N_EQS() == solver_data->get_N_EQS());
    assert(U.size() == config.get_N_CELLS_INT() && R.size() == config.get_N_CELLS_INT());

    /*--------------------------------------------------------------------
     U_n+1 = U_n + dt /Omega * R(U_n)
    --------------------------------------------------------------------*/

    evaluate_flux_balance(config, U);
    for_all(U, i, n_eq)
        U(i, n_eq) += dt / cells.get_cell_volume(i) * R(i, n_eq);
}

void Solver::TVD_RK3(const Config &config)
{
    const Scalar dt = config.get_delta_time();
    VecField &U = solver_data->get_solution();
    VecField &U_old = solver_data->get_solution_old();
    VecField &R = solver_data->get_flux_balance();
    const auto &cells = FV_grid.get_cells();
    Index i, n_eq;

    assert(U.get_N_EQS() == solver_data->get_N_EQS() && U_old.get_N_EQS() == solver_data->get_N_EQS() && R.get_N_EQS() == solver_data->get_N_EQS());
    assert(U.size() == config.get_N_CELLS_INT() && U_old.size() == config.get_N_CELLS_INT());
    assert(R.size() == config.get_N_CELLS_INT());

    /*--------------------------------------------------------------------
    U_1 = U_n + dt / Omega * R(U_n)
    U_2 = 3/4 * U_n + 1/4 *U_1 + 1/4 * dt / Omega * R(U_1)
    U_n+1 = 1/3 * U_n + 2/3 * U_2 + 2/3 * dt / Omega * R(U_2)
    --------------------------------------------------------------------*/

    evaluate_flux_balance(config, U);
    for_all(U, i, n_eq)
        U(i, n_eq) += dt / cells.get_cell_volume(i) * R(i, n_eq);

    evaluate_flux_balance(config, U);
    for_all(U, i, n_eq)
        U(i, n_eq) = 3.0 / 4.0 * U_old(i, n_eq) + 1.0 / 4.0 * U(i, n_eq) + 1.0 / 4.0 * dt / cells.get_cell_volume(i) * R(i, n_eq);

    evaluate_flux_balance(config, U);
    for_all(U, i, n_eq)
        U(i, n_eq) = 1.0 / 3.0 * U_old(i, n_eq) + 2.0 / 3.0 * U(i, n_eq) + 2.0 / 3.0 * dt / cells.get_cell_volume(i) * R(i, n_eq);
}

void Solver::communicate_primvars()
{
    part_comm.clear();
    part_comm.pack_field(solver_data->get_primvars());
    part_comm.communicate_ghost_fields();
    part_comm.unpack_field(solver_data->get_primvars());
}
void Solver::communicate_primvars_and_gradient()
{
    part_comm.clear();
    part_comm.pack_field(solver_data->get_primvars());
    part_comm.pack_field(solver_data->get_primvars_gradient());
    part_comm.communicate_ghost_fields();
    part_comm.unpack_field(solver_data->get_primvars());
    part_comm.unpack_field(solver_data->get_primvars_gradient());
}
void Solver::communicate_max_and_min()
{
    part_comm.clear();
    part_comm.pack_field(solver_data->get_primvars_max());
    part_comm.pack_field(solver_data->get_primvars_min());
    part_comm.communicate_ghost_fields();
    part_comm.unpack_field(solver_data->get_primvars_max());
    part_comm.unpack_field(solver_data->get_primvars_min());
}
void Solver::communicate_limiter()
{
    part_comm.clear();
    part_comm.pack_field(solver_data->get_primvars_limiter());
    part_comm.communicate_ghost_fields();
    part_comm.unpack_field(solver_data->get_primvars_limiter());
}

EulerSolver::EulerSolver(const Config &config,
                         const geometry::PrimalGrid &primal_grid,
                         const geometry::FV_Grid &FV_grid,
                         PartitionComm &part_comm)
    : Solver(config, primal_grid, FV_grid, part_comm)
{
    solver_data = make_unique<EulerSolverData>(config);
    validity_checker = make_unique<EulerValidityChecker>(config);

    calc_Delta_S(config);

    ShortIndex max_scalars_per_cell =
        (solver_data->get_n_vecfields_sendrecv_max() + solver_data->get_n_gradfields_sendrecv_max() * N_DIM) * solver_data->get_N_EQS();
    part_comm.set_size(max_scalars_per_cell);
}

void EulerSolver::evaluate_inviscid_fluxes(const Config &config)
{
    VecField &flux_balance = solver_data->get_flux_balance();

    // const VecField& primvars = solver_data->get_primvars();
    // const GradField& primvars_grad = solver_data->get_primvars_gradient();
    // const VecField& primvars_limiter = solver_data->get_primvars_limiter();

    Index N_INTERIOR_FACES = config.get_N_FACES_INT();
    SpatialOrder spatial_order = config.get_spatial_order();

    Index i, j;
    const auto &faces = FV_grid.get_faces();
    const auto &patches = FV_grid.get_patches_boundary();

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
        i = faces.get_cell_i(ij);
        j = faces.get_cell_j(ij);
        const Vec3 &S_ij = faces.get_face_normal(ij);
        const Vec3 &r_im = faces.get_centroid_to_face_i(ij);
        const Vec3 &r_jm = faces.get_centroid_to_face_j(ij);

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
    This is handled PatchBoundary-wise*/

    for (Index i_PatchBoundary{0}; i_PatchBoundary < patches.size(); i_PatchBoundary++)
    {
        const auto &PatchBoundary = patches[i_PatchBoundary];
        auto &boundary_condition = BC_container[i_PatchBoundary];
        // BoundaryCondition::BC_function BC_func = BoundaryCondition::get_BC_function(PatchBoundary.boundary_type);

        for (Index ij{PatchBoundary.FIRST_FACE}; ij < PatchBoundary.FIRST_FACE + PatchBoundary.N_FACES; ij++)
        {

            Index i_domain = faces.get_cell_i(ij);
            const Vec3 &S_ij = faces.get_face_normal(ij);
            const Vec3 &r_im = faces.get_centroid_to_face_i(ij);

            calc_reconstructed_value(i_domain, V_L, r_im, spatial_order);

            boundary_condition->calc_ghost_val(V_L, V_R, S_ij);

            EulerEqs::prim_to_cons(V_L, U_L);
            EulerEqs::prim_to_cons(V_R, U_R);

            numerical_flux_func(U_L, U_R, S_ij, Flux_inv);

            assert(validity_checker->valid_boundary_flux(Flux_inv.data(), PatchBoundary.boundary_type));

            flux_balance.get_variable<EulerVec>(i_domain) -= Flux_inv;
        }
    }
    /*Finally interprocessor communication*/
}

void EulerSolver::calc_timestep(Config &config)
{
    // --------------------------------------------------------------------
    // Implementing Method 2 in "Time Step on Unstructured Grids" from Blazek
    // --------------------------------------------------------------------

    /*Delta S only needs updating when the grid is moved*/
    if (config.get_grid_motion())
        calc_Delta_S(config);

    const EulerSolverData &euler_data = dynamic_cast<const EulerSolverData &>(*solver_data);
    const vector<Vec3> &Delta_S = euler_data.get_Delta_S();

    const Scalar CFL = config.get_CFL();
    auto cells = FV_grid.get_cells();

    const auto &primvars = solver_data->get_primvars();

    Scalar c, volume, spec_rad_x, spec_rad_y, spec_rad_z;

    Scalar delta_time = std::numeric_limits<Scalar>::max(); // Large number

    for (Index i{0}; i < config.get_N_CELLS_INT(); i++)
    {
        const EulerVecMap V = primvars.get_variable<EulerVec>(i);
        c = EulerEqs::sound_speed_primitive(V);

        volume = cells.get_cell_volume(i);

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

    const auto &faces = FV_grid.get_faces();
    EulerSolverData &euler_data = dynamic_cast<EulerSolverData &>(*solver_data);

    vector<Vec3> &Delta_S = euler_data.get_Delta_S();
    std::fill(Delta_S.begin(), Delta_S.end(), Vec3::Zero());
    Index i, j;
    Vec3 tmp;

    for (Index ij{0}; ij < config.get_N_FACES_INT(); ij++)
    {
        i = faces.get_cell_i(ij);
        j = faces.get_cell_j(ij);
        const Vec3 &S_ij = faces.get_face_normal(ij);

        tmp = 0.5 * S_ij.cwiseAbs();

        Delta_S[i] += tmp;
        Delta_S[j] += tmp;
    }
}

void EulerSolver::set_constant_ghost_values(const Config &config)
{
    const auto &patches = FV_grid.get_patches_boundary();
    const auto &faces = FV_grid.get_faces();
    VecField &primvars = solver_data->get_primvars();
    Index i_domain, j_ghost;

    // for (const auto &PatchBoundary : patches)
    for (Index i_PatchBoundary{0}; i_PatchBoundary < patches.size(); i_PatchBoundary++)
    {
        const auto &PatchBoundary = patches[i_PatchBoundary];
        auto &boundary_condition = BC_container[i_PatchBoundary];

        for (Index ij{PatchBoundary.FIRST_FACE}; ij < PatchBoundary.FIRST_FACE + PatchBoundary.N_FACES; ij++)
        {
            i_domain = faces.get_cell_i(ij);
            j_ghost = faces.get_cell_j(ij);
            assert(i_domain < config.get_N_CELLS_INT() && j_ghost >= config.get_N_CELLS_INT());

            const Vec3 &S_ij = faces.get_face_normal(ij);
            // primvars.map_to_variable<EulerVec>(j) =
            //     BoundaryCondition::calc_ghost_val<BC_type, EulerVec>(primvars.map_to_variable<EulerVec>(i), S_ij);

            // Test if this works
            EulerVecMap V_i_domain = primvars.get_variable<EulerVec>(i_domain);
            EulerVecMap V_j_ghost = primvars.get_variable<EulerVec>(j_ghost);

            boundary_condition->calc_ghost_val(V_i_domain, V_j_ghost, S_ij);
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
        gradient::calc_green_gauss_gradient<N_EQS_EULER>(config, FV_grid, primvars, primvars_grad);
        break;
    default:
        assert(false); // no others are yet implemented
    }
}

void EulerSolver::evaluate_limiter(const Config &config)
{
    assert(config.get_limiter() != Limiter::NONE);
    const VecField &primvars = solver_data->get_primvars();
    const GradField &primvars_grad = solver_data->get_primvars_gradient();
    VecField &primvars_limiter = solver_data->get_primvars_limiter();
    VecField &primvars_max = solver_data->get_primvars_max();
    VecField &primvars_min = solver_data->get_primvars_min();

    switch (config.get_limiter())
    {
    case Limiter::Barth:
        communicate_max_and_min();
        reconstruction::calc_max_and_min_values<N_EQS_EULER>(config,
                                                             FV_grid,
                                                             primvars,
                                                             primvars_max,
                                                             primvars_min);

        reconstruction::calc_barth_limiter<N_EQS_EULER>(config,
                                                        FV_grid,
                                                        primvars,
                                                        primvars_grad,
                                                        primvars_max,
                                                        primvars_min,
                                                        primvars_limiter);
        break;
    default:
        FAIL_MSG("Selected limiter not implemented\n");
    }
    communicate_limiter();
}
