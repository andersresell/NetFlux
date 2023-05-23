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

    void evaluate_flux_balance(Config& config);

    double calc_timestep(Config& config) override;

    SolverType get_solver_type() const override {return SolverType::Euler;}
private:

    void evaluate_gradient(const Config& config);
    void evaluate_limiter(const Config& config);

    template<InviscidFluxScheme Scheme>
    void evaluate_inviscid_flux_balance(const Config& config);

    static void linear_reconstruct_face(Index i,
                                        const VecField& primvars,
                                        const VecField& primvars_grad,
                                        const VecField& primvars_limiter,
                                        const Vec3& r_im,
                                        EulerVec& V_face
                                        );
    static void constant_reconstruct_face(Index i,
                                          const VecField& primvars,
                                          const EulerVec* V_face);
};




template<InviscidFluxScheme Scheme>
void EulerSolver::evaluate_inviscid_flux_balance(const Config& config){
    
    VecField& flux_balance = solver_data->get_flux_balance();
    const VecField& primvars = solver_data->get_primvars();
    const GradField& primvars_grad = solver_data->get_primvars_gradient();
    const VecField& primvars_limiter = solver_data->get_primvars_limiter();

    Index N_INTERIOR_FACES = config.get_N_INTERIOR_FACES();
    Index N_TOTAL_FACES = config.get_N_TOTAL_FACES();
    SpatialOrder spatial_order = config.get_spatial_order();

    Index i,j;
    const auto& cells = grid.get_cells();
    const auto& faces = grid.get_faces();
    const auto& patches = grid.get_patches();

    //InvFluxFunction inv_flux_func = NumericalFlux::get_inviscid_flux_function(config);
    EulerVec Flux_inv; //inviscid numerical flux
    Vec3 normal;
    double face_area;

    EulerVec V_L, V_R, U_L, U_R;

    /*First interior cells*/
    for (Index ij{0}; ij<N_INTERIOR_FACES; ij++){
        i = faces[ij].i;
        j = faces[ij].j;
        const Vec3& S_ij = faces[ij].S_ij;
        const Vec3& r_im = faces[ij].r_im;
        const Vec3& r_jm = faces[ij].r_jm;
        
        if (spatial_order == SpatialOrder::Second){
            
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
        
        NumericalFlux::inviscid_flux<Scheme>(U_L, U_R, S_ij, Flux_inv);

        flux_balance[i] -= Flux_inv;
        flux_balance[j] += Flux_inv;
    }

    /*Then boundaries. Here ghost cells has to be assigned based on the boundary conditions.
    This is handled patch-wise*/
    for (const auto& patch : patches){

        for (Index ij{patch.FIRST_FACE}; ij<patch.FIRST_FACE+patch.N_FACES; ij++){
            if (spatial_order == SpatialOrder::Second){
            
            auto V_i = static_cast<const EulerVec&>(primvars[i]);
            auto V_i_grad = static_cast<const EulerGrad&>(primvars_grad[i]);
            auto limiter_i = static_cast<const EulerVec&>(primvars_limiter[i]);
            Reconstruction::calc_reconstructed_value<N_EQS_EULER>(V_i, V_i_grad, limiter_i, r_im, V_L);

        } else{ //First order spatial
            V_L = static_cast<const EulerVec&>(primvars[i]);
            V_R = static_cast<const EulerVec&>(primvars[j]);
        }
        }
    }


}


void linear_reconstruct_face(Index i,
                            const VecField& primvars,
                            const VecField& primvars_grad,
                            const VecField& primvars_limiter,
                            const Vec3& r_im,
                            EulerVec* V_face
                            ){
    auto V_i = static_cast<const EulerVec&>(primvars[i]);
    auto V_i_grad = static_cast<const EulerGrad&>(primvars_grad[i]);
    auto limiter_i = static_cast<const EulerVec&>(primvars_limiter[i]);
    Reconstruction::calc_reconstructed_value<N_EQS_EULER>(V_i, V_i_grad, limiter_i, r_im, V_L);
}
void constant_reconstruct_face(Index i,
                                const VecField& primvars,
                                EulerVec* V_face){
    V_face = primvars.cast<EulerVec>(i);
    *V_face = static_cast<EulerVec>(primvars[i]);
}



class NS_Solver : public EulerSolver{

};