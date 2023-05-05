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
    void zero_residual();

};













class NS_Solver : public EulerSolver{

};