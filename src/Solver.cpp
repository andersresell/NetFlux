#include "../include/Solver.hpp"

using namespace geom;

Solver::Solver(const geom::Grid& grid)
    : grid{grid}
{

}

EulerSolver::EulerSolver(const Config& config, const geom::Grid& grid) : Solver(grid)
{
    solver_data = make_unique<EulerSolverData>(config);

}


void EulerSolver::step(Config& config){
    
    assert(config.get_time_integration_type() == TimeIntegrationType::Explicit); //Remove if implementing implicit

    double dt = config.get_delta_time();

}

void EulerSolver::evaluate_residual(Config& config){
    VecField& residual = solver_data->get_residual();
    residual.set_zero();

    if (config.get_spatial_order() == SpatialOrder::Second){
        evaluate_gradient(config);

        if (config.get_limiter() != Limiter::NONE)
            evaluate_limiter(config);
    }

    switch (config.get_inv_flux_scheme()){
        case InviscidFluxScheme::Rusanov:
            inviscid_flux_balance<InviscidFluxScheme::Rusanov>(config);
            break;
        default:
            FAIL_MSG("Illegal numerical inviscid flux scheme\n");
    }

}

double EulerSolver::calc_timestep(Config& config){
    // --------------------------------------------------------------------
    // Implementing Method 2 in "Time Step on Unstructured Grids" in Blazek 
    // --------------------------------------------------------------------
    
    const double CFL = config.get_CFL();
    auto cells = grid.get_cells();
    auto faces = grid.get_faces();

    const auto& U = solution->cell_values;
    
    double rho, u, v, w, c, volume, spec_rad_x, spec_rad_y, spec_rad_z;

    double delta_time = std::numeric_limits<double>::max(); //Large number

    for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
        //solution[i].sound_speed();  
        rho = U[i][0];
        u = U[i][1]/rho;
        v = U[i][2]/rho;
        w = U[i][3]/rho;
        c = EulerField::sound_speed(U[i]);
        volume = cells[i].cell_volume;

        double Delta_S_x{0}, Delta_S_y{0}, Delta_S_z{0}; //These are projections of the control volume in each spatial direction
        
        const auto& neigbour_faces = grid.get_surrounding_faces(i);

        for (Index ij : neigbour_faces){
            const Vec3& S_ij = faces[ij].S_ij;
            Delta_S_x += abs(S_ij.x());
            Delta_S_y += abs(S_ij.y());
            Delta_S_z += abs(S_ij.z());
        }
        Delta_S_x *= 0.5;
        Delta_S_y *= 0.5;
        Delta_S_z *= 0.5;

        spec_rad_x = (abs(u) + c) * Delta_S_x;
        spec_rad_y = (abs(v) + c) * Delta_S_y;
        spec_rad_z = (abs(w) + c) * Delta_S_z;
                
        delta_time = std::min(delta_time, CFL* volume / (spec_rad_x + spec_rad_y + spec_rad_z));
    }
    config.set_delta_time(delta_time);
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
