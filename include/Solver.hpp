#pragma once

#include "../include/Config.hpp"
#include "../include/Grid.hpp"
#include "../include/Utilities.hpp"
#include "../include/Numerics.hpp"
#include "../include/SolverData.hpp"


class Solver{
protected:
    unique_ptr<SolverData> solver_data;
    const geom::Grid& grid;

public:

    Solver(const geom::Grid& grid);

    virtual void step(Config& config) = 0;
    virtual double calc_timestep(Config& config) = 0;
    virtual SolverType get_solver_type() const = 0;

    const SolverData& get_solver_data() const {return *solver_data;} 

    

};


class EulerSolver : public Solver{   

public:
    EulerSolver(const Config& config, const geom::Grid& grid);

    void step(Config& config) override;

    void evaluate_residual(Config& config);


    double calc_timestep(Config& config) override;

    SolverType get_solver_type() const override {return SolverType::Euler;}
private:

    void evaluate_gradient(const Config& config);
    void evaluate_limiter(const Config& config);

    template<InviscidFluxScheme Scheme>
    void inviscid_flux_balance(const Config& config);

};




template<InviscidFluxScheme Scheme>
void EulerSolver::inviscid_flux_balance(const Config& config){

    VecField& primvars = solver_data->get_primvars();
    GradField& primvars_grad = solver_data->get_primvars_gradient();
    VecField& primvars_limiter = solver_data->get_primvars_limiter();

    Index i,j;
    const auto& cells = grid.get_cells();
    const auto& faces = grid.get_faces();
    InvFluxFunction inv_flux_func = NumericalFlux::get_inviscid_flux_function(config);
    EulerVec Flux_inv_ij; //inviscid numerical flux
    Vec3 normal;
    double face_area;

    EulerVec V_L, V_R, U_L, U_R;

    for (Index ij{0}; ij<config.get_N_FACES(); ij++){
        i = faces[ij].i;
        j = faces[ij].j;
        const Vec3& S_ij = faces[ij].S_ij;
        const Vec3& r_im = faces[ij].r_im;
        const Vec3& r_jm = faces[ij].r_jm;
        
        if (config.get_spatial_order() == SpatialOrder::Second){
            
            auto V_i = static_cast<const EulerVec&>(primvars[i]);
            auto V_i_grad = static_cast<const EulerGrad&>(primvars_grad[i]);
            auto limiter_i = static_cast<const EulerVec&>(primvars_limiter[i]);
            Reconstruction::calc_reconstructed_value<N_EQS_EULER>(V_i, V_i_grad, limiter_i, r_im, V_L);


            auto V_j = static_cast<const EulerVec&>(primvars[j]);
            auto V_j_grad = static_cast<const EulerGrad&>(primvars_grad[j]);
            auto limiter_j = static_cast<const EulerVec&>(primvars_limiter[j]);
            Reconstruction::calc_reconstructed_value<N_EQS_EULER>(V_j, V_j_grad, limiter_j, r_jm, V_R);

        } else{ //First order spatial
            V_L = static_cast<const EulerVec&>(primvars[i]);
            V_R = static_cast<const EulerVec&>(primvars[j]);
        }


        U_L = EulerEqs::prim_to_cons(V_L);
        U_R = EulerEqs::prim_to_cons(U_R);
        
        NumericalFlux::inviscid_flux<Scheme>(U_L, U_R, S_ij, Flux_inv_ij);
    }
}









class NS_Solver : public EulerSolver{

};