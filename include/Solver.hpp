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

    void evaluate_flux_balance(Config& config); //move to Solver class later

    double calc_timestep(Config& config) override;

    SolverType get_solver_type() const override {return SolverType::Euler;}
private:

    void set_constant_ghost_values(const Config& config);
    void evaluate_gradient(const Config& config);
    void evaluate_limiter(const Config& config);

    void inviscid_flux_balance(const Config& config);

    void calc_reconstructed_value(Index i, 
                                  EulerVecMap& V_L,  
                                  const Vec3& r_im, 
                                  SpatialOrder spatial_order);

    void calc_ghost_value(const EulerVecMap& U_L, EulerVecMap& U_R, const Vec3& S_ij, BoundaryType boundary_type);
};







class NS_Solver : public EulerSolver{

};