#pragma once

#include "../include/Config.hpp"
#include "../include/Grid.hpp"
#include "../include/Utilities.hpp"
#include "../include/Numerics.hpp"
#include "../include/SolverData.hpp"

class Solver
{
protected:
    unique_ptr<SolverData> solver_data;
    const geom::Grid &grid;
    unique_ptr<ValidityChecker> validity_checker;

    Vector<unique_ptr<BoundaryCondition>> BC_container;

public:
    Solver(const geom::Grid &grid, const Config &config);

    void step(const Config &config);

    virtual void calc_timestep(Config &config) = 0;

    virtual SolverType get_solver_type() const = 0;

    const SolverData &get_solver_data() const { return *solver_data; }

private:
    void create_BC_container(const Config &config);

    void evaluate_flux_balance(const Config &config, const VecField &cons_vars);

    void update_cell_averages(VecField &U, double dt);

    virtual void evaluate_inviscid_fluxes(const Config &config) { return; };

    virtual void evaluate_viscous_fluxes(const Config &config) { return; };

    virtual void set_constant_ghost_values(const Config &config) = 0;

    void explicit_euler(const Config &config);

    void TVD_RK3(const Config &config);

    virtual void evaluate_gradient(const Config &config) = 0;

    virtual void evaluate_limiter(const Config &config) = 0;
};

class EulerSolver : public Solver
{

public:
    EulerSolver(const Config &config, const geom::Grid &grid);

    void calc_timestep(Config &config) override;

    SolverType get_solver_type() const override { return SolverType::Euler; }

private:
    void evaluate_inviscid_fluxes(const Config &config) final;

    void set_constant_ghost_values(const Config &config) final;

    void evaluate_gradient(const Config &config) final;

    void evaluate_limiter(const Config &config) final;

    void calc_reconstructed_value(Index i,
                                  EulerVecMap &V_L,
                                  const Vec3 &r_im,
                                  SpatialOrder spatial_order);

    /*Delta S is used to compute time step following the 2nd method in Blazek.
    Only needs recalculating when the grid is updated*/
    void calc_Delta_S(const Config &config);
};

inline void EulerSolver::calc_reconstructed_value(Index i,
                                                  EulerVecMap &V_L,
                                                  const Vec3 &r_im,
                                                  SpatialOrder spatial_order)
{

    const VecField &primvars = solver_data->get_primvars();
    const GradField &primvars_grad = solver_data->get_primvars_gradient();
    const VecField &primvars_limiter = solver_data->get_primvars_limiter();

    const EulerVecMap V_i = primvars.get_variable<EulerVec>(i);

    if (spatial_order == SpatialOrder::Second)
    {
        const EulerGradMap V_i_grad = primvars_grad.get_variable<EulerGrad>(i);
        const EulerVecMap limiter_i = primvars_limiter.get_variable<EulerVec>(i);

        Reconstruction::calc_limited_reconstruction(V_i, V_i_grad, limiter_i, r_im, V_L);
    }
    else
    { // Spatial order == First
        V_L = V_i;
    }
}

class NS_Solver : public EulerSolver
{
};