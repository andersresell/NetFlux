#pragma once

#include "../include/Config.hpp"
#include "../include/Grid.hpp"
#include "../include/Utilities.hpp"
#include "../include/Numerics.hpp"


class BaseSolver{

};

template<ShortIndex N_EQS>
class Solver : public BaseSolver{

protected:

    FlowField<N_EQS> U, //solution 
                     U_old, //old solution
                     R; //residual

    const geom::Grid& grid;
};



class EulerSolver : public Solver<N_EQS_EULER>{
    using EulerVar = flow::EulerVar;
    

public:
    const Vector<EulerVar>& get_solution() const {return U;}

    EulerSolver(const Config& config, const geom::Grid& grid);

    void step(Config& config);

    void evaluate_residual(Config& config);


    double calc_timestep(Config& config);

private:
    void zero_residual();

};













class NS_Solver : public EulerSolver{

};