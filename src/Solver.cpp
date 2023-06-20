#include "../include/Solver.hpp"

using namespace geom;

Solver::Solver(const geom::Grid& grid)
    : grid{grid}
{

}

void Solver::step(const Config& config){
    
    assert(config.get_time_integration_type() == TimeIntegrationType::Explicit); //Remove if implementing implicit

    switch (config.get_time_scheme()){
        case TimeScheme::ExplicitEuler:
            explicit_euler(config);
            break;
        case TimeScheme::TVD_RK3:
            TVD_RKD(config);
            break;
        default:
            FAIL_MSG("Error, illegal time integration scheme used\n");
    }

    solver_data->get_solution_old() = solver_data->get_solution();

}

void Solver::evaluate_flux_balance(const Config& config, const VecField& cons_vars){
    
    solver_data->get_flux_balance().set_zero();
    
    solver_data->set_primvars(cons_vars, config);

    set_constant_ghost_values(config);

    if (config.get_spatial_order() == SpatialOrder::Second){
        evaluate_gradient(config);

        if (config.get_limiter() != Limiter::NONE)
            evaluate_limiter(config);
    }

    evaluate_inviscid_fluxes(config);

    evaluate_viscous_fluxes(config);

}

void Solver::explicit_euler(const Config& config){
    double dt = config.get_delta_time();
    VecField& U = solver_data->get_solution();
    VecField& R = solver_data->get_flux_balance();
    const ShortIndex N_EQS = solver_data->get_N_EQS();
    assert(U.size() == R.size() && U.get_N_EQS() == N_EQS && R.get_N_EQS() == N_EQS);
    
    //Modify evaluate_flux_balance and following functions, so that they are evaluated based on a certain field U or 
    //perhaps privars V is sufficient to determine everything. 
    evaluate_flux_balance(config, U);
    
    for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
        for (ShortIndex j{0}; j<N_EQS; j++){
            U(i,j) += dt * R(i,j);
        }
    }

}

void Solver::TVD_RKD(const Config& config){
    // double dt = config.get_delta_time();
    // VecField& U = solver_data->get_solution();
    // VecField& R = solver_data->get_flux_balance();
    // const ShortIndex N_EQS = solver_data->get_N_EQS();
    // assert(U.size() == R.size() && U.get_N_EQS() == N_EQS && R.get_N_EQS() == N_EQS);
    
    // evaluate_flux_balance(config);
    
    // for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
    //     for (ShortIndex j{0}; j<N_EQS; j++){
    //         U(i,j) += dt * R(i,j);
    //     }
    // }
}



EulerSolver::EulerSolver(const Config& config, const geom::Grid& grid) : Solver(grid)
{
    solver_data = make_unique<EulerSolverData>(config);


    calc_Delta_S(config);
}





void EulerSolver::evaluate_inviscid_fluxes(const Config& config){
    
    VecField& flux_balance = solver_data->get_flux_balance();
    //const VecField& primvars = solver_data->get_primvars();
    //const GradField& primvars_grad = solver_data->get_primvars_gradient();
    //const VecField& primvars_limiter = solver_data->get_primvars_limiter();

    Index N_INTERIOR_FACES = config.get_N_INTERIOR_FACES();
    SpatialOrder spatial_order = config.get_spatial_order();

    Index i,j;
    const auto& faces = grid.get_faces();
    const auto& patches = grid.get_patches();

    //InvFluxFunction inv_flux_func = NumericalFlux::get_inviscid_flux_function(config);


    EulerSolverData& euler_data = dynamic_cast<EulerSolverData&>(*solver_data);
    EulerVecMap U_L = euler_data.get_U_L_map();
    EulerVecMap U_R = euler_data.get_U_R_map();
    EulerVecMap V_L = euler_data.get_V_L_map();
    EulerVecMap V_R = euler_data.get_V_R_map();
    EulerVecMap Flux_inv = euler_data.get_Flux_inv_map();

    NumericalFlux::InvFluxFunction inviscid_flux = NumericalFlux::get_inviscid_flux_function(config);

    /*First interior cells*/
    for (Index ij{0}; ij<N_INTERIOR_FACES; ij++){
        i = faces[ij].i;
        j = faces[ij].j;
        const Vec3& S_ij = faces[ij].S_ij;
        const Vec3& r_im = faces[ij].r_im;
        const Vec3& r_jm = faces[ij].r_jm;
        
        
        calc_reconstructed_value(i, V_L, r_im, spatial_order);
        calc_reconstructed_value(j, V_R, r_jm, spatial_order);    

        EulerEqs::prim_to_cons(V_L, U_L);
        EulerEqs::prim_to_cons(V_R, U_R);
        
        inviscid_flux(U_L, U_R, S_ij, Flux_inv);

        flux_balance.get_variable<EulerVec>(i) -= Flux_inv;
        flux_balance.get_variable<EulerVec>(j) += Flux_inv;
    }

    /*Then boundaries. Here ghost cells has to be assigned based on the boundary conditions.
    This is handled patch-wise*/

    for (const auto& patch : patches){
        
        BoundaryCondition::BC_function BC_func = BoundaryCondition::get_BC_function(patch.boundary_type);

        for (Index ij{patch.FIRST_FACE}; ij<patch.FIRST_FACE+patch.N_FACES; ij++){

            i = faces[ij].i;
            j = faces[ij].j;
            const Vec3& S_ij = faces[ij].S_ij;
            const Vec3& r_im = faces[ij].r_im;

            calc_reconstructed_value(i, V_L, r_im, spatial_order);
            
            //calc_ghost_value(V_L, V_R, S_ij, patch.boundary_type);
            
            BC_func(V_L, V_R, S_ij);

            EulerEqs::prim_to_cons(V_L, U_L);
            EulerEqs::prim_to_cons(V_R, U_R);
            
            inviscid_flux(U_L, U_R, S_ij, Flux_inv);

            flux_balance.get_variable<EulerVec>(i) -= Flux_inv;
      
        }
    }
}

void EulerSolver::calc_reconstructed_value(Index i, 
                                            EulerVecMap& V_L,  
                                            const Vec3& r_im, 
                                            SpatialOrder spatial_order){

    const VecField& primvars = solver_data->get_primvars();
    const GradField& primvars_grad = solver_data->get_primvars_gradient();
    const VecField& primvars_limiter = solver_data->get_primvars_limiter();
    
    const EulerVecMap V_i = primvars.get_variable<EulerVec>(i);

    if (spatial_order == SpatialOrder::Second){
        const EulerGradMap V_i_grad = primvars_grad.get_variable<EulerGrad>(i);
        const EulerVecMap limiter_i = primvars_limiter.get_variable<EulerVec>(i);

        Reconstruction::calc_limited_reconstruction(V_i, V_i_grad, limiter_i, r_im, V_L);

    } else{ //Spatial order == First
        V_L = V_i;
    }
}




// double EulerSolver::calc_timestep(Config& config){
//     // --------------------------------------------------------------------
//     // Implementing Method 2 in "Time Step on Unstructured Grids" in Blazek 
//     // --------------------------------------------------------------------
    
//     const double CFL = config.get_CFL();
//     auto cells = grid.get_cells();
//     auto faces = grid.get_faces();

//     const auto& U = solution->cell_values;
    
//     double rho, u, v, w, c, volume, spec_rad_x, spec_rad_y, spec_rad_z;

//     double delta_time = std::numeric_limits<double>::max(); //Large number

//     for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
//         //solution[i].sound_speed();  
//         rho = U[i][0];
//         u = U[i][1]/rho;
//         v = U[i][2]/rho;
//         w = U[i][3]/rho;
//         c = EulerField::sound_speed(U[i]);
//         volume = cells[i].cell_volume;

//         double Delta_S_x{0}, Delta_S_y{0}, Delta_S_z{0}; //These are projections of the control volume in each spatial direction
        
//         const auto& neigbour_faces = grid.get_surrounding_faces(i);

//         for (Index ij : neigbour_faces){
//             const Vec3& S_ij = faces[ij].S_ij;
//             Delta_S_x += abs(S_ij.x());
//             Delta_S_y += abs(S_ij.y());
//             Delta_S_z += abs(S_ij.z());
//         }
//         Delta_S_x *= 0.5;
//         Delta_S_y *= 0.5;
//         Delta_S_z *= 0.5;

//         spec_rad_x = (abs(u) + c) * Delta_S_x;
//         spec_rad_y = (abs(v) + c) * Delta_S_y;
//         spec_rad_z = (abs(w) + c) * Delta_S_z;
                
//         delta_time = std::min(delta_time, CFL* volume / (spec_rad_x + spec_rad_y + spec_rad_z));
//     }
//     config.set_delta_time(delta_time);
// }

void EulerSolver::calc_timestep(Config& config){
    // --------------------------------------------------------------------
    // Implementing Method 2 in "Time Step on Unstructured Grids" in Blazek 
    // --------------------------------------------------------------------
    
    /*Delta S only needs updating when the grid is moved*/
    if (config.get_grid_motion()) 
        calc_Delta_S(config);

    const EulerSolverData& euler_data = dynamic_cast<const EulerSolverData&>(*solver_data);
    const Vector<Vec3>& Delta_S = euler_data.get_Delta_S();

    const double CFL = config.get_CFL();
    auto cells = grid.get_cells();

    const auto& primvars = solver_data->get_primvars();
    
    double c, volume, spec_rad_x, spec_rad_y, spec_rad_z;

    double delta_time = std::numeric_limits<double>::max(); //Large number

    for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
        const EulerVecMap V = primvars.get_variable<EulerVec>(i);
        c = EulerEqs::sound_speed_primitive(V);

        volume = cells[i].cell_volume;

        spec_rad_x = (abs(V[1]) + c) * Delta_S[i].x();
        spec_rad_y = (abs(V[2]) + c) * Delta_S[i].y();
        spec_rad_z = (abs(V[3]) + c) * Delta_S[i].z();

        delta_time = std::min(delta_time, CFL* volume / (spec_rad_x + spec_rad_y + spec_rad_z));
    }
    config.set_delta_time(delta_time);
}

void EulerSolver::calc_Delta_S(const Config& config){
    const auto& faces = grid.get_faces();
    EulerSolverData& euler_data = dynamic_cast<EulerSolverData&>(*solver_data);
    
    Vector<Vec3>& Delta_S = euler_data.get_Delta_S();
    std::fill(Delta_S.begin(), Delta_S.end(), Vec3::Zero());
    Index i,j;
    Vec3 tmp;

    for (Index ij{0}; ij<config.get_N_INTERIOR_FACES(); ij++){
        i = faces[ij].i;
        j = faces[ij].j;

        tmp = 0.5*Delta_S[ij].cwiseAbs();
        
        Delta_S[i] += tmp;
        Delta_S[j] += tmp;
    }
}

void EulerSolver::set_constant_ghost_values(const Config& config){

    const auto& patches = grid.get_patches();
    const auto& faces = grid.get_faces();
    VecField& primvars = solver_data->get_primvars();
    Index i,j;

    EulerVecMap V_i{nullptr};
    EulerVecMap V_j{nullptr};

    for (const auto& patch : patches){
    
        BoundaryCondition::BC_function bc_func = BoundaryCondition::get_BC_function(patch.boundary_type);

        for (Index ij{patch.FIRST_FACE}; ij<patch.FIRST_FACE+patch.N_FACES; ij++){
            i = faces[ij].i;
            j = faces[ij].j;
            const Vec3& S_ij = faces[ij].S_ij;
            // primvars.map_to_variable<EulerVec>(j) = 
            //     BoundaryCondition::calc_ghost_val<BC_type, EulerVec>(primvars.map_to_variable<EulerVec>(i), S_ij);

            //Test if this works
            V_i = primvars.get_variable<EulerVec>(i);
            V_j = primvars.get_variable<EulerVec>(j);
             
            bc_func(V_i, V_j, S_ij);
            
        }
    }
}

void EulerSolver::evaluate_gradient(const Config& config){
    const VecField& primvars= solver_data->get_primvars();
    GradField& primvars_grad = solver_data->get_primvars_gradient();
    
    switch(config.get_grad_scheme()){
        case GradientScheme::GreenGauss:
            Gradient::calc_green_gauss_gradient<N_EQS_EULER>(config, grid, primvars, primvars_grad);
            break;
        default:
            FAIL_MSG("Selected gradient scheme not implemented\n");
    }
}   


void EulerSolver::evaluate_limiter(const Config& config){

    const VecField& primvars= solver_data->get_primvars();
    const GradField& primvars_grad = solver_data->get_primvars_gradient();
    VecField& primvars_limiter = solver_data->get_primvars_limiter();
    VecField& primvars_max = solver_data->get_primvars_max();
    VecField& primvars_min = solver_data->get_primvars_min();


    switch(config.get_limiter()){
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

