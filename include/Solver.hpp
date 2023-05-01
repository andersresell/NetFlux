#pragma once

#include "../include/Config.hpp"
#include "../include/Grid.hpp"
#include "../include/Utilities.hpp"
#include "../include/Numerics.hpp"

/*template<typename FlowVariable>
class Solver{
protected:
    virtual FlowVariable calc_flux(const FlowVariable& U) const = 0;
};*/



class EulerSolver {
    using EulerVar = flow::EulerVar;
    
    Vector<EulerVar> U, //solution 
                     U_old, //old solution
                     R; //residual

    const geom::Grid& grid;

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