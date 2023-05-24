#pragma once

#include "Utilities.hpp"
#include "Grid.hpp"
#include "SolverData.hpp"

//DEPRECATED FUNCTION POINTER IMPLEMENTATION
/*
using InvFluxFunction = std::function<void(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux)>;

class NumericalFlux{

public:

    static void Rusanov(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux);

    static void HLLC(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux);

public:
    static InvFluxFunction get_inviscid_flux_function(const Config& config);

};*/


namespace NumericalFlux{
    /*Templated riemann solvers*/
    using namespace geom;

    template<InviscidFluxScheme Scheme>
    inline void inviscid_flux(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux);

    /*Rusanov implementation*/
    template<>
    inline void inviscid_flux<InviscidFluxScheme::Rusanov>(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux){
        
        Vec3 normal = S_ij.normalized();
        double area = S_ij.norm();
        double spec_rad_L = EulerEqs::conv_spec_rad(U_L, normal);
        double spec_rad_R = EulerEqs::conv_spec_rad(U_R, normal);


        Flux = area * 0.5*(EulerEqs::inviscid_flux(U_R, normal) + EulerEqs::inviscid_flux(U_L, normal) 
            - std::max(spec_rad_R, spec_rad_L) * (U_R - U_L));
            
    }
    
}




namespace Gradient{
    using namespace geom;

    /*Implementing the "compact" gradient in the cell center, from chapter 9.2 in Moukalled et. al. No orthogonal correction for now*/

    template<ShortIndex N_EQS>
    inline void calc_green_gauss_gradient( const Config& config,
                        const Grid& grid,
                        const VecField& vec_field,
                        GradField& grad_field){
                            
        assert(grad_field.cols() == N_DIM && grad_field.rows() == vec_field.rows());
        assert(N_EQS == vec_field.get_N_EQS());
    
        using FieldVec = StaticContainer1D<double, N_EQS>;
        using FieldGrad = StaticContainer2D<double, N_EQS, N_DIM>;


        const auto& faces = grid.get_faces();
        const auto& cells = grid.get_cells();

        const Index N_FACES = config.get_N_FACES();
        const Index N_CELLS = config.get_N_INTERIOR_CELLS();
        Index i,j;
        
        FieldVec U_Face;
        FieldGrad tmp;

        grad_field.set_zero();

        for (Index ij{0}; ij<N_FACES; ij++){
            const Face& face = faces[ij];

            const Cell& cell_i = cells[face.i];
            const Cell& cell_j = cells[face.j];
            i = cell_i.i;
            j = cell_i.j;

            

            //U_face = 0.5 * (vec_field.get_eigen_matrix(i) + vec_field.get_eigen_matrix(j));
            U_face = 0.5 * ((FieldVec)vec_field[i] + (FieldVec)vec_field[j]);

            //Simple average for now, might improve later with distance weighting and orthogonal correctors later
            container::Vec3_mult((FieldVec)vec_field[ij], face.S_ij, tmp);
            grad_field[i] += tmp / cell_i.cell_volume;
            if (j < N_CELLS) //only calculate gradient for interior cells?
                grad_field[j] -= tmp / cell_j.cell_volume;      
        }
        
    }
}

namespace Reconstruction{
    using namespace geom;



    template<ShortIndex N_EQS>
    inline void calc_reconstructed_value(
                                   const FlowVec<N_EQS>& V_c,
                                   const FlowGrad<N_EQS>& V_c_grad,
                                   const FlowVec<N_EQS>& limiter_c,
                                   const Vec3& r_cf, 
                                   FlowVec<N_EQS>& V_f
                                   ){
        
        for (ShortIndex k{0}; k<N_EQS; k++){
            double Delta_V{0};
            for (ShortIndex iDim{0}; iDim<N_DIM; iDim++)
                Delta_V += V_c_grad(k, iDim) * r_cf[iDim];

            V_f[k] = V_c[k] + limiter_c[k] * Delta_V[k];
        }
    }

    /*Implementing the Barth limiter procedure in Blazek*/
    template<ShortIndex N_EQS> 
    inline void calc_barth_limiter( const Config& config,
                        const Grid& grid,
                        const VecField& sol_field,
                        const GradField& sol_grad,
                        const VecField& max_field,
                        const VecField& min_field,
                        VecField& limiter){
        
        assert(N_EQS == sol_field.get_N_EQS() && N_EQS == sol_grad.get_N_EQS() && N_EQS == max_field.get_N_EQS() && 
            N_EQS == min_field.get_N_EQS() && N_EQS == limiter.get_N_EQS());
        assert(sol_field.length() == sol_grad.length() == max_field.length() == min_field.length() == limiter.length());


        const auto& faces = grid.get_faces();
        const auto& cells = grid.get_cells();

        const Index N_FACES = config.get_N_FACES();
        const Index N_CELLS = config.get_N_INTERIOR_CELLS();
        Index i,j;
        Eigen::VectorXd Delta_2;
        gradient.resize(N_EQS);
        double EPS = std::numeric_limits<double>::epsilon(); 

        limiter = DBL_MAX;

        for (Index ij{0}; ij<N_FACES; ij++){
            const geom::Face& face = faces[ij];

            i = cells[face.i].i;
            j = cells[face.j].j;

            const Eigen::MatrixX3d& gradient_i = sol_grad.get_eigen_matrix(i);
            Delta_2 = gradient_i * face.r_im;
            assert(Delta_2.rows() == N_EQS);

            for (ShortIndex k{0}; k<N_EQS; k++){
                //to avoid division by zero    
                Delta_2[k] = sign(Delta_2[k])*(abs(Delta_2[k]) + EPS);

                if (Delta_2[k] > 0.0)
                    limiter(i,k) = min(limiter(i,k), min(1.0, (max_field(i,k) - sol_field(i,k))/Delta_2));
                else if (Delta_2 < 0.0)
                    limiter(i,k) = min(limiter(i,k), min(1.0, (min_field(i,k) - sol_field(i,k))/Delta_2));
                else
                    limiter(i,k) = min(limiter(i,k), 1.0);
            }

            if (j < N_CELLS){
                const Eigen::MatrixX3d& gradient_j = sol_grad.get_eigen_matrix(j);
                Delta_2 = gradient_j * face.r_jm;

                for (ShortIndex k{0}; k<N_EQS; k++){
                    //to avoid division by zero    
                    Delta_2[k] = sign(Delta_2[k])*(abs(Delta_2[k]) + EPS);

                    if (Delta_2[k] > 0.0)
                        limiter(j,k) = min(limiter(j,k), min(1.0, (max_field(j,k) - sol_field(j,k))/Delta_2));
                    else if (Delta_2 < 0.0)
                        limiter(j,k) = min(limiter(j,k), min(1.0, (min_field(j,k) - sol_field(j,k))/Delta_2));
                    else
                        limiter(j,k) = min(limiter(j,k), 1.0);
                }

            }

        } 
    }





    /*Used in calculation of limiters*/
    template<ShortIndex N_EQS> 
    void calc_max_and_min_values( const Config& config,
                        const Grid& grid,
                        const VecField& sol_field,
                        VecField& max_field,
                        VecField& min_field){

        assert(N_EQS == sol_field.get_N_EQS() && N_EQS == max_field.get_N_EQS() && N_EQS == min_field.get_N_EQS());
        assert(sol_field.length() == max_field.length() == min_field.length());

        /*Setting the max and min fields to either a very large or a very small value*/
        max_field = -DBL_MAX;
        min_field = DBL_MAX;
    

        const auto& faces = grid.get_faces();
        const auto& cells = grid.get_cells();

        const Index N_FACES = config.get_N_FACES();
        const Index N_CELLS = config.get_N_INTERIOR_CELLS();
        Index i,j;
        
        /*the formulas to be computed are 
        U_max = max(U_i, max_j(U_j)) and U_min = min(U_i, min_j(U_j))
        Looping over all edges instead of cells due to efficiency*/

        for (Index ij{0}; ij<N_FACES; ij++){
            const Face& face = faces[ij];

            i = cells[face.i].i;
            j = cells[face.j].j;

            for (ShortIndex k{0}; k<N_EQS; k++){
                max_field(i,k) = std::max(max_field(i,k), sol_field(i,k));
                max_field(i,k) = std::max(max_field(i,k), sol_field(j,k));
                if (j < N_CELLS) //Only for interior cells
                    max_field(j,k) = std::max(max_field(j,k), sol_field(i,k));

                min_field(i,k) = std::min(min_field(i,k), sol_field(i,k));
                min_field(i,k) = std::min(min_field(i,k), sol_field(j,k));
                if (j < N_CELLS) //Only for interior cells
                    min_field(j,k) = std::min(min_field(j,k), sol_field(i,k));
            }

        }    
    }
}
