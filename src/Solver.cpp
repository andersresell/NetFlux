#include "../include/Solver.hpp"

using namespace flow;
using namespace geom;

EulerSolver::EulerSolver(const Config& config, const geom::Grid& grid)
    : grid{grid}
{
    //Specify intial val rather
    U.resize(config.get_N_INTERIOR_CELLS());
    U_old.resize(config.get_N_INTERIOR_CELLS());
    R.resize(config.get_N_INTERIOR_CELLS());
    
}


void EulerSolver::step(Config& config){
    
    assert(config.get_time_integration_type() == TimeIntegrationType::Explicit); //Remove if implementing implicit

    double dt = config.get_delta_time();

}

void EulerSolver::evaluate_residual(Config& config){
    zero_residual();

    Index i,j;
    const auto& cells = grid.get_cells();
    const auto& faces = grid.get_faces();
    InvFluxFunction inv_flux_func = NumericalFlux::get_inviscid_flux_function(config);
    EulerVar F_inv_ij; //inviscid numerical flux

    for (Index ij{0}; ij<config.get_N_FACES(); ij++){
        i = faces[ij].i;
        j = faces[ij].j;
        const Vec3& S_ij = faces[ij].S_ij;
        
        //Calculate gradients etc. 

        const auto& U_L; //recinstructed state at cell i side;
        const auto& U_R; //reconstructed state at cell j side;
        
        
        inv_flux_func(U_L, U_R, F_inv_ij);
    }

}

double EulerSolver::calc_timestep(Config& config){
    // --------------------------------------------------------------------
    // Implementing Method 2 in "Time Step on Unstructured Grids" in Blazek 
    // --------------------------------------------------------------------
    
    
    assert(config.get_governing_eq() == GoverningEq::Euler); //edit this when adding NavierStokes

    const double CFL = config.get_CFL();
    auto cells = grid.get_cells();
    auto faces = grid.get_faces();
    
    double rho, u, v, w, c, volume, spec_rad_x, spec_rad_y, spec_rad_z;

    double delta_time = std::numeric_limits<double>::max(); //Large number

    for (Index i{0}; i<config.get_N_INTERIOR_CELLS(); i++){
        //solution[i].sound_speed();  
        rho = U[i][0];
        u = U[i][1]/rho;
        v = U[i][2]/rho;
        w = U[i][3]/rho;
        c = EulerVar::sound_speed(U[i]);
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

void EulerSolver::zero_residual(){
    for (Index i{0}; i<R.size(); i++)
        for (Index j{0}; j<N_EQS_EULER; j++)
            R[i][j] = 0.0;
}