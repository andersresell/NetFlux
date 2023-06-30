#pragma once
#include "Includes.hpp"
#include "containers/DynamicContainer.hpp"
#include "containers/StaticContainer.hpp"
#include "Config.hpp"

struct VecField final : public DynamicContainer2D<Scalar>
{
    VecField(Index N_CELLS, ShortIndex N_EQS) : DynamicContainer2D(N_CELLS, N_EQS) {}

    ShortIndex get_N_EQS() const { return rows(); }

    using DynamicContainer2D<Scalar>::operator=;
};

struct GradField final : public DynamicContainer3D<Scalar, N_DIM>
{
    GradField(Index N_CELLS, ShortIndex N_EQS) : DynamicContainer3D(N_CELLS, N_EQS) {}

    ShortIndex get_N_EQS() const { return rows(); }

    using DynamicContainer3D<Scalar, N_DIM>::operator=;
};

class SolverData
{
protected:
    unique_ptr<VecField> solution,
        solution_old,
        primvars,
        flux_balance;

    unique_ptr<GradField> primvars_gradient;

    unique_ptr<VecField> primvars_limiter,
        primvars_max,
        primvars_min;

    // SolverData() = default;

    SolverData(const Config &config, ShortIndex n_eqs);

public:
    virtual ShortIndex get_N_EQS() const = 0;

    VecField &get_solution() { return *solution; }
    const VecField &get_solution() const { return *solution; }

    VecField &get_solution_old() { return *solution_old; }
    const VecField &get_solution_old() const { return *solution_old; }

    VecField &get_flux_balance() { return *flux_balance; }
    const VecField &get_flux_balance() const { return *flux_balance; }

    VecField &get_primvars() { return *primvars; }
    const VecField &get_primvars() const { return *primvars; }

    GradField &get_primvars_gradient() { return *primvars_gradient; }
    const GradField &get_primvars_gradient() const { return *primvars_gradient; }

    VecField &get_primvars_limiter() { return *primvars_limiter; }
    const VecField &get_primvars_limiter() const { return *primvars_limiter; }

    VecField &get_primvars_max() { return *primvars_max; }
    const VecField &get_primvars_max() const { return *primvars_max; }

    VecField &get_primvars_min() { return *primvars_min; }
    const VecField &get_primvars_min() const { return *primvars_min; }

    virtual SolverType get_solver_type() const = 0;

    virtual void set_primvars(const VecField &cons_vars, const Config &config) = 0;
};

template <ShortIndex N_EQS>
using FlowVec = Eigen::Vector<Scalar, N_EQS>;
template <ShortIndex N_EQS>
using FlowGrad = Eigen::Matrix<Scalar, N_EQS, N_DIM>;

// using EulerVec = StaticContainer1D<Scalar, N_EQS_EULER>;
// using EulerGrad = StaticContainer2D<Scalar, N_EQS_EULER, N_DIM>;
using EulerVec = FlowVec<N_EQS_EULER>;
using EulerGrad = FlowGrad<N_EQS_EULER>;

using EulerVecMap = Eigen::Map<EulerVec>;
using EulerGradMap = Eigen::Map<EulerGrad>;

class EulerSolverData : public SolverData
{

    EulerVec U_L, U_R, V_L, V_R;
    EulerVec Flux_inv;

    Vector<Vec3> Delta_S; // Used in time step computation

public:
    EulerSolverData(const Config &config);

    ShortIndex get_N_EQS() const final { return N_EQS_EULER; };

    virtual SolverType get_solver_type() const { return SolverType::Euler; }

    EulerVecMap get_U_L_map() { return EulerVecMap(U_L.data()); }
    EulerVecMap get_U_R_map() { return EulerVecMap(U_R.data()); }
    EulerVecMap get_V_L_map() { return EulerVecMap(V_L.data()); }
    EulerVecMap get_V_R_map() { return EulerVecMap(V_R.data()); }
    EulerVecMap get_Flux_inv_map() { return EulerVecMap(Flux_inv.data()); }

    EulerVecMap get_empty_map() const { return EulerVecMap(nullptr); }

    Vector<Vec3> &get_Delta_S() { return Delta_S; }
    const Vector<Vec3> &get_Delta_S() const { return Delta_S; }

    void set_primvars(const VecField &cons_vars, const Config &config) final;

    void set_freestream_values(const Config &config);
};

/*Discontinuing StaticContainer, using Eigen instead*/

// template<ShortIndex N_EQS>
// using FlowVec = StaticContainer1D<Scalar, N_EQS>;
// template<ShortIndex N_EQS>
// using FlowGrad = StaticContainer2D<Scalar, N_EQS, N_DIM>;

namespace EulerEqs
{

    constexpr Scalar GAMMA{standard_air::gamma};
    constexpr Scalar GAMMA_MINUS_ONE{GAMMA - 1};
    constexpr Scalar GAMMA_MINUS_ONE_INV{1 / GAMMA_MINUS_ONE};

    /*Templating the various functions since the EulerVecType can be either EulerVec or EulerVecMap*/

    template <typename EulerVecType>
    inline void prim_to_cons(const EulerVecType &V, EulerVecType &U);

    template <typename EulerVecType>
    inline void cons_to_prim(const EulerVecType &U, EulerVecType &V);

    template <typename EulerVecType>
    inline Scalar pressure(const EulerVecType &U);

    template <typename EulerVecType>
    inline Scalar sound_speed_conservative(const EulerVecType &U);

    template <typename EulerVecType>
    inline Scalar sound_speed_primitive(const EulerVecType &V);

    template <typename EulerVecType>
    inline Scalar projected_velocity(const EulerVecType &U, const Vec3 &normal);

    template <typename EulerVecType>
    inline Scalar conv_spectral_radii(const EulerVecType &U, const Vec3 &normal);

    template <typename EulerVecType>
    inline Scalar conv_spectral_radii(const EulerVecType &U, const Vec3 &normal);

    template <typename EulerVecType>
    inline EulerVec inviscid_flux(const EulerVecType &U, const Vec3 &normal);

    template <typename EulerVecType>
    inline void prim_to_cons(const EulerVecType &V, EulerVecType &U)
    {
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER and EulerVecType::ColsAtCompileTime == 1);
        U[0] = V[0];
        U[1] = V[0] * V[1];
        U[2] = V[0] * V[2],
        U[3] = V[0] * V[3],
        U[4] = GAMMA_MINUS_ONE_INV * V[4] + 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2] + V[3] * V[3]);

        assert(!U.hasNaN() && U.allFinite());
        assert(U[4] > 0.0);
    }

    template <typename EulerVecType>
    inline void cons_to_prim(const EulerVecType &U, EulerVecType &V)
    {

        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER and EulerVecType::ColsAtCompileTime == 1);
        V[0] = U[0],
        V[1] = U[1] / U[0],
        V[2] = U[2] / U[0],
        V[3] = U[3] / U[0],
        V[4] = pressure(U);

        assert(!V.hasNaN() && V.allFinite());
        assert(V[0] > 0.0 && V[4] > 0.0);
    }

    template <typename EulerVecType>
    inline Scalar pressure(const EulerVecType &U)
    {
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        Scalar p = GAMMA_MINUS_ONE * (U[4] - 0.5 / U[0] * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]));
        assert(num_is_valid_and_pos(p));
        return p;
    }

    template <typename EulerVecType>
    inline Scalar sound_speed_conservative(const EulerVecType &U)
    {
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        Scalar c = sqrt(GAMMA * pressure(U) / U[0]);
        assert(num_is_valid_and_pos(c));
        return c;
    }

    template <typename EulerVecType>
    inline Scalar sound_speed_primitive(const EulerVecType &V)
    {
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        Scalar c = sqrt(GAMMA * V[4] / V[0]);
        assert(num_is_valid_and_pos(c));
        return c;
    }

    template <typename EulerVecType>
    inline Scalar projected_velocity(const EulerVecType &U, const Vec3 &normal)
    {
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        Scalar vel_n = (U[1] * normal.x() + U[2] * normal.y() + U[3] * normal.z()) / U[0];
        assert(num_is_valid(vel_n));
        return vel_n;
    }

    template <typename EulerVecType>
    inline Scalar conv_spectral_radii(const EulerVecType &U, const Vec3 &normal)
    {
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        Scalar lambda = abs(projected_velocity(U, normal)) + sound_speed_conservative(U);
        assert(num_is_valid_and_pos(lambda));
        return lambda;
    }

    template <typename EulerVecType>
    inline EulerVec inviscid_flux(const EulerVecType &U, const Vec3 &normal)
    {
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);

        Scalar rho = U[0];
        Scalar p = pressure(U);
        Scalar V_normal = (U[1] * normal.x() + U[2] * normal.y() + U[3] * normal.z()) / rho;

        EulerVec F = {V_normal * rho,
                      V_normal * U[1] + p * normal.x(),
                      V_normal * U[2] + p * normal.y(),
                      V_normal * U[3] + p * normal.z(),
                      V_normal * (U[4] + p)};
        assert(!F.hasNaN() && F.allFinite());
        return F;
    }

}

/*Implements various mechanisms for checking if the various containers (solutions, fluxes etc) contain physical solutions*/
class ValidityChecker
{
protected:
    const Config &config;
    /*Checks for inf or nan values*/
    Index check_field_validity(const VecField &field, Index first, Index last) const;

    /*Special check for primvars (for instance, for Euler eqs. density and pressure must be positive)*/
    virtual Index check_primvars(const VecField &V, Index first, Index last) const = 0;

    virtual Index check_consvars(const VecField &U, Index first, Index last) const = 0;

    virtual string get_solver_name() const = 0;

    void write_debug_info(const VecField &U) const;

public:
    void check_flux_balance_validity(const Config &config, const VecField &flux_balance) const;

    bool valid_primvars_interior(const VecField &V) const;

    bool valid_consvars_interior(const VecField &U) const;

    bool valid_primvars_ghost(const VecField &V) const;

    bool valid_consvars_ghost(const VecField &U) const;

    bool valid_flux_balance(const VecField &R) const;

    ValidityChecker(const Config &config) : config{config} {}
};

class EulerValidityChecker : public ValidityChecker
{

    Index check_primvars(const VecField &V, Index first, Index last) const final;

    Index check_consvars(const VecField &U, Index first, Index last) const final;

    string get_solver_name() const override { return "Euler"; }

public:
    EulerValidityChecker(const Config &config) : ValidityChecker(config) {}
};
