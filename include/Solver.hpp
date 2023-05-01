#pragma once

#include "../include/Config.hpp"
#include "../include/Grid.hpp"
#include "../include/Utilities.hpp"

/*template<typename FlowVariable>
class Solver{
protected:
    virtual FlowVariable calc_flux(const FlowVariable& U) const = 0;
};*/



class EulerSolver {
    using EulerVar = flow::EulerVar;
    
    Vector<EulerVar> solution, solution_old;

    const geom::Grid& grid;

public:
    const Vector<EulerVar>& get_solution() const {return solution;}

    EulerSolver(const Config& config, const geom::Grid& grid);

    void step(Config& config);




    double calc_timestep(Config& config);

};













class NS_Solver : public EulerSolver{

};