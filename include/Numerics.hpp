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
    template<class VecType, class GradType>
    void calc_green_gauss_gradient(const Config& config,
                     const Grid& grid,
                     const Vector<VecType>& vec_field,
                     Vector<GradType>& grad_field);



    template<class VecType, class GradType>
    inline void calc_green_gauss_gradient( const Config& config,
                        const Grid& grid,
                        const Vector<VecType>& vec_field,
                        Vector<GradType>& grad_field){
        static_assert(grad_field.rows() == N_DIM && grad_field.columns() == vec_field.size());

        const auto& faces = grid.get_faces();
        const auto& cells = grid.get_cells();

        const Index N_FACES = config.get_N_FACES();
        const Index N_CELLS = config.get_N_INTERIOR_CELLS();
        Index i,j;
        
        GradType tmp;
        VecType face_value;

        for (Index i_cell = 0; i_cell<grad_field.size(); i_cell++)
            grad_field[i_cell] = 0; //zeroing grad field
            
        for (Index ij{0}; ij<N_FACES){
            const Face& face = faces[ij];

            const Cell& cell_i = cells[face.i];
            const Cell& cell_j = cells[face.j];
            i = cell_i.i;
            j = cell_i.j;

            face_value = 0.5 * (vec_field[i] + vec_field[j]);  

            //Simple average for now, might improve later with distance weighting and orthogonal correctors later
            container::Vec3_mult(vec_field[ij], face.S_ij, tmp);
            grad_field[i] += tmp / cell_i.cell_volume;
            if (j < N_CELLS) //only calculate gradient for interior cells?
                grad_field[j] -= tmp / cell_j.cell_volume;      
        }
        
    }
}
