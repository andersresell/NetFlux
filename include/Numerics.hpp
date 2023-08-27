#pragma once

#include "Utilities.hpp"
#include "geometry/FV_Grid.hpp"
#include "SolverData.hpp"

class NumericalFlux
{
public:
    using InvFluxFunction = std::function<void(const EulerVecMap &U_L, const EulerVecMap &U_R, const Vec3 &S_ij, EulerVecMap &Flux)>;

    static InvFluxFunction get_inviscid_flux_function(InviscidFluxScheme inv_flux_scheme);

private:
    static void rusanov(const EulerVecMap &U_L, const EulerVecMap &U_R, const Vec3 &S_ij, EulerVecMap &Flux);

    static void HLLC(const EulerVecMap &U_L, const EulerVecMap &U_R, const Vec3 &S_ij, EulerVecMap &Flux);

    static inline const map<InviscidFluxScheme, InvFluxFunction> inv_flux_functions = {{InviscidFluxScheme::Rusanov, rusanov},
                                                                                       {InviscidFluxScheme::HLLC, HLLC}};
};

// class BoundaryCondition
// {
// public:
//     using BC_function = std::function<void(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij)>;

//     static BC_function get_BC_function(BoundaryType boundary_type);

// private:
//     static void no_slip_wall(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij);

//     static void slip_wall(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij);

//     static void farfield(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij);

//     static inline const map<BoundaryType, BC_function> BC_functions = {{BoundaryType::NoSlipWall, no_slip_wall},
//                                                                        {BoundaryType::SlipWall, slip_wall},
//                                                                        {BoundaryType::FarField, farfield}};
// };

class BoundaryCondition
{

public:
    virtual void calc_ghost_val(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij) = 0;
};

class BC_NoSlipWall : public BoundaryCondition
{
public:
    void calc_ghost_val(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij) final;
};

class BC_SlipWall : public BoundaryCondition
{
    Vec3 normal;
    Scalar vel_normal;

public:
    void calc_ghost_val(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij) final;
};

class BC_FarField : public BoundaryCondition
{
    Scalar c_domain, c_fs, c_boundary;
    Scalar entropy_fs, entropy_domain, entropy_boundary;
    Scalar density_fs, density_boundary;
    Vec3 vel_domain, vel_fs, vel_boundary;
    Scalar pressure_fs, pressure_boundary;
    Vec3 normal;
    Scalar vel_n_domain, vel_n_fs, vel_n_boundary;
    Scalar Riemann_plus, Riemann_minus;

public:
    BC_FarField(const Config &config);
    void calc_ghost_val(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij) final;
};

namespace gradient
{
    using namespace geometry;

    /*Implementing the "compact" gradient in the cell center, from chapter 9.2 in Moukalled et. al. No orthogonal correction for now*/

    template <ShortIndex N_EQS>
    inline void calc_green_gauss_gradient(const Config &config,
                                          const FV_Grid &FV_grid,
                                          const VecField &vec_field,
                                          GradField &grad_field)
    {

        assert(grad_field.cols() == N_DIM && grad_field.rows() == N_EQS && vec_field.rows() == N_EQS);
        assert(grad_field.size() == config.get_N_CELLS_DOMAIN() && vec_field.size() == config.get_N_CELLS_TOT());

        using FieldVec = Eigen::Vector<Scalar, N_EQS>;
        using FieldGrad = Eigen::Matrix<Scalar, N_EQS, N_DIM>;

        const Faces &faces = FV_grid.get_faces();
        const Cells &cells = FV_grid.get_cells();

        const Index N_FACES = config.get_N_FACES_TOT();
        const Index N_CELLS = config.get_N_CELLS_INT();

        Index i, j;

        FieldVec U_face;
        FieldGrad tmp;

        grad_field.set_zero();

        for (Index ij{0}; ij < N_FACES; ij++)
        {
            i = faces.get_cell_i(ij);
            j = faces.get_cell_j(ij);
            const Vec3 &S_ij = faces.get_face_normal(ij);

            // Simple average for now, might improve with distance weighting and orthogonal correctors later
            U_face = 0.5 * (vec_field.get_variable<FieldVec>(i) + vec_field.get_variable<FieldVec>(j));

            tmp = U_face * S_ij.transpose(); // DOES THIS MAKE SENSE?? (ij)

            grad_field.get_variable<FieldGrad>(i) += tmp / cells.get_cell_volume(i);

            if (j < N_CELLS) // only calculate gradient for interior cells?
                grad_field.get_variable<FieldGrad>(j) -= tmp / cells.get_cell_volume(j);
        }
    }
}

namespace reconstruction
{
    using namespace geometry;

    template <typename VecMapType, typename GradMapType, typename VecType>
    inline void calc_limited_reconstruction(
        const VecMapType &V_center,
        const GradMapType &V_center_grad,
        const VecMapType &limiter_center,
        const Vec3 &r_center2face,
        VecType &V_face)
    {

        V_face = V_center + limiter_center.cwiseProduct(V_center_grad * r_center2face);
    }

    /*Implementing the Barth limiter procedure in Blazek*/
    template <ShortIndex N_EQS>
    inline void calc_barth_limiter(const Config &config,
                                   const FV_Grid &FV_grid,
                                   const VecField &sol_field,
                                   const GradField &sol_grad,
                                   const VecField &max_field,
                                   const VecField &min_field,
                                   VecField &limiter)
    {

        assert(N_EQS == sol_field.get_N_EQS() && N_EQS == sol_grad.get_N_EQS() && N_EQS == max_field.get_N_EQS() &&
               N_EQS == min_field.get_N_EQS() && N_EQS == limiter.get_N_EQS());

        const Index N_CELLS_INTERIOR = config.get_N_CELLS_INT();
        const Index N_TOTAL_FACES = config.get_N_FACES_TOT();

        assert(sol_field.size() == config.get_N_CELLS_TOT() && sol_grad.size() == config.get_N_CELLS_DOMAIN() && max_field.size() == config.get_N_CELLS_DOMAIN() &&
               min_field.size() == config.get_N_CELLS_DOMAIN() && limiter.size() == config.get_N_CELLS_DOMAIN());

        using FieldVec = Eigen::Vector<Scalar, N_EQS>;
        using FieldGrad = Eigen::Matrix<Scalar, N_EQS, N_DIM>;
        using FieldGradMap = Eigen::Map<FieldGrad>;

        constexpr Scalar EPS = std::numeric_limits<Scalar>::epsilon();

        const Faces &faces = FV_grid.get_faces();

        FieldVec Delta_2;
        FieldGrad gradient;
        Scalar face_contribution;

        limiter = DBL_MAX;

        for (Index ij{0}; ij < N_TOTAL_FACES; ij++)
        {
            Index i = faces.get_cell_i(ij);
            Index j = faces.get_cell_j(ij);
            const Vec3 &r_im = faces.get_centroid_to_face_i(ij);
            const Vec3 &r_jm = faces.get_centroid_to_face_j(ij);
            const FieldGradMap gradient_i = sol_grad.get_variable<FieldGrad>(i);
            Delta_2 = gradient_i * r_im;

            for (ShortIndex k{0}; k < N_EQS; k++)
            {
                Delta_2[k] = sign(Delta_2[k]) * (abs(Delta_2[k]) + EPS);

                if (Delta_2[k] > 0.0)
                    face_contribution = min(1.0, (max_field(i, k) - sol_field(i, k)) / Delta_2[k]);
                else if (Delta_2[k] < 0.0)
                    face_contribution = min(1.0, (min_field(i, k) - sol_field(i, k)) / Delta_2[k]);
                else
                    face_contribution = 1.0;

                // Branchless version
                // face_contribution = (Delta_2[k] > 0.0) * min(1.0, (max_field(i, k) - sol_field(i, k)) / Delta_2[k]) +
                //                     (Delta_2[k] < 0.0) * min(1.0, (min_field(i, k) - sol_field(i, k)) / Delta_2[k]) +
                //                     (Delta_2[k] < 0.0) * 1.0;

                limiter(i, k) = min(limiter(i, k), face_contribution);
                assert(limiter(i, k) >= 0.0 && limiter(i, k) <= 1.0);
            }

            if (j < N_CELLS_INTERIOR)
            {
                const FieldGradMap &gradient_j = sol_grad.get_variable<FieldGrad>(j);
                Delta_2 = gradient_j * r_jm;

                for (ShortIndex k{0}; k < N_EQS; k++)
                {
                    // to avoid division by zero
                    Delta_2[k] = sign(Delta_2[k]) * (abs(Delta_2[k]) + EPS);

                    if (Delta_2[k] > 0.0)
                        face_contribution = min(1.0, (max_field(j, k) - sol_field(j, k)) / Delta_2[k]);
                    else if (Delta_2[k] < 0.0)
                        face_contribution = min(1.0, (min_field(j, k) - sol_field(j, k)) / Delta_2[k]);
                    else
                        face_contribution = 1.0;

                    // Branchless version
                    // face_contribution = (Delta_2[k] > 0.0) * min(1.0, (max_field(j, k) - sol_field(j, k)) / Delta_2[k]) +
                    //                     (Delta_2[k] < 0.0) * min(1.0, (min_field(j, k) - sol_field(j, k)) / Delta_2[k]) +
                    //                     (Delta_2[k] == 0.0) * 1.0;

                    limiter(j, k) = min(limiter(j, k), face_contribution);
                    assert(limiter(j, k) >= 0.0 && limiter(j, k) <= 1.0);
                }
            }
        }

#ifndef NDEBUG
        /*Checking that the values lay between 0 and 1*/
        constexpr Scalar TOL = 1e-8;
        for (Index i{0}; i < N_CELLS_INTERIOR; i++)
        {
            for (Index j = 0; j < limiter.rows(); j++)
            {
                assert(num_is_valid(limiter(i, j)));
                assert(limiter(i, j) > -TOL && limiter(i, j) < 1.0 + TOL);
            }
        }
#endif
    }

    /*Used in calculation of limiters*/
    template <ShortIndex N_EQS>
    void calc_max_and_min_values(const Config &config,
                                 const FV_Grid &FV_grid,
                                 const VecField &sol_field,
                                 VecField &max_field,
                                 VecField &min_field)
    {

        const Index N_TOTAL_FACES = config.get_N_FACES_TOT();
        const Index N_INTERIOR_CELLS = config.get_N_CELLS_INT();
        assert(N_EQS == sol_field.get_N_EQS() && N_EQS == max_field.get_N_EQS() && N_EQS == min_field.get_N_EQS());
        assert(sol_field.size() == config.get_N_CELLS_TOT());
        assert(max_field.size() == config.get_N_CELLS_DOMAIN());
        assert(min_field.size() == config.get_N_CELLS_DOMAIN());

        /*Setting the max and min fields to either a very large or a very small value*/
        max_field = -DBL_MAX;
        min_field = DBL_MAX;

        const auto &faces = FV_grid.get_faces();

        /*the formulas to be computed are (Taken from Blazek)
        U_max = max(U_i, max_j(U_j)) and U_min = min(U_i, min_j(U_j))
        Looping over all edges instead of cells due to efficiency*/

        for (Index ij{0}; ij < N_TOTAL_FACES; ij++)
        {
            Index i = faces.get_cell_i(ij);
            Index j = faces.get_cell_j(ij);

            for (ShortIndex k{0}; k < N_EQS; k++)
            {
                max_field(i, k) = max(max_field(i, k), max(sol_field(i, k), sol_field(j, k)));
                if (j < N_INTERIOR_CELLS) // Only for interior cells
                    max_field(j, k) = max(max_field(j, k), max(sol_field(i, k), sol_field(j, k)));

                min_field(i, k) = min(min_field(i, k), min(sol_field(i, k), sol_field(j, k)));
                if (j < N_INTERIOR_CELLS) // Only for interior cells
                    min_field(j, k) = min(min_field(j, k), min(sol_field(i, k), sol_field(j, k)));
            }
        }
    }
}
