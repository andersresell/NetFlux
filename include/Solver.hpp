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

    void step(const Config& config);

    virtual double calc_timestep(Config& config) = 0;

protected:

    const SolverData& get_solver_data() const {return *solver_data;} 

private:   

    virtual void evaluate_flux_balance(const Config& config, const VecField& conservative) = 0;

    virtual SolverType get_solver_type() const = 0;

    void explicit_euler(const Config& config);

    void TVD_RKD(const Config& config);


};


class EulerSolver : public Solver{   

public:
    EulerSolver(const Config& config, const geom::Grid& grid);

    double calc_timestep(Config& config) override;

    SolverType get_solver_type() const override {return SolverType::Euler;}

private:
    void evaluate_flux_balance(const Config& config, const VecField& cons_vars) override; 

    void set_constant_ghost_values(const Config& config);
    void evaluate_gradient(const Config& config);
    void evaluate_limiter(const Config& config);

    void inviscid_flux_balance(const Config& config);

    void calc_reconstructed_value(Index i, 
                                  EulerVecMap& V_L,  
                                  const Vec3& r_im, 
                                  SpatialOrder spatial_order);
    
    /*Delta S is used to compute time step following the 2nd method in Blazek. 
    Only needs recalculating when the grid is updated*/
    void calc_Delta_S(const Config& config);


};







class NS_Solver : public EulerSolver{

};