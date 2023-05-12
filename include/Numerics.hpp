#pragma once

#include "Utilities.hpp"
#include "Grid.hpp"
#include "SolverData.hpp"

using InvFluxFunction = std::function<void(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux)>;

class NumericalFlux{

public:

    static void Rusanov(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux);

    static void HLLC(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux);

public:
    static InvFluxFunction get_inviscid_flux_function(const Config& config);

};




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
    
        using FieldVec = StaticStaticContainer1D<double, N_EQS>;
        using FieldGrad = StaticContainer2D<double, N_EQS, N_DIM>;


        const auto& faces = grid.get_faces();
        const auto& cells = grid.get_cells();

        const Index N_FACES = config.get_N_FACES();
        const Index N_CELLS = config.get_N_INTERIOR_CELLS();
        Index i,j;
        
        //ShortIndex N_EQS = vec_field.get_N_EQS();
        
        FieldVec U_Face;
        FieldGrad tmp;

        //Mat3X tmp; //  U * S_ij
        //tmp.resize(N_EQS);
        //VecX U_face; //the value of the field U interpolated at the face;
        //U_face.resize(N_EQS)

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
