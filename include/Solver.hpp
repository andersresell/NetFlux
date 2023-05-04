#pragma once

#include "../include/Config.hpp"
#include "../include/Grid.hpp"
#include "../include/Utilities.hpp"
#include "../include/Numerics.hpp"
#include "../include/FlowVariable.hpp"


class BaseSolver{
public:
    virtual void step(Config& config) = 0;
    virtual double calc_timestep(Config& config);
};

template<ShortIndex N_EQS>
class Solver : public BaseSolver{
protected:

    unique_pointer<FlowField<N_EQS>> solution, 
                     solution_old,
                     residual;

    
    const geom::Grid& grid;
public:
    const FlowField& get_solution() const {return *solution;} 

    Solver(const geom::Grid& grid);

    

};



class EulerSolver : public Solver<N_EQS_EULER>{
    
    

public:
    EulerSolver(const Config& config, const geom::Grid& grid);

    void step(Config& config) override;

    void evaluate_residual(Config& config);


    double calc_timestep(Config& config) override;


private:
    void zero_residual();

};













class NS_Solver : public EulerSolver{

};