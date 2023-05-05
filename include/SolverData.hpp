#pragma once
#include "includes.hpp"
#include "Container.hpp"
#include "Config.hpp"




/*template<Index N_EQS>
struct FlowVar : public Container1D<double, N_EQS> {
    FlowVar() = default;
    FlowVar(std::initializer_list<double> init) : Container1D(init) {} 
};*/

constexpr double GAMMA{1.4};
constexpr double GAMMA_MINUS_ONE{1-GAMMA};
constexpr double GAMMA_MINUS_ONE_INV{1/GAMMA_MINUS_ONE};

/*Variables of the Euler / NS equations, can be both conservative, primitive, fluxes, etc.*/
// struct EulerVar : public Container1D<double, N_EQS_EULER>{
//     EulerVar() =  default;
    
//     //EulerVar(std::initializer_list<double> init) : Container1D(init) {}
//     template<typename... Args>
//     EulerVar(Args&&... args) : Container1D(std::forward<Args>(args)...) {}
    
//     static EulerVar prim_to_cons(const EulerVar& V);
//     static EulerVar cons_to_prim(const EulerVar& U);

//     static double pressure(const EulerVar& U);
//     static double sound_speed(const EulerVar& U);

//     static EulerVar inviscid_flux_x(const EulerVar& U);
//     static EulerVar inviscid_flux_y(const EulerVar& U);
//     static EulerVar inviscid_flux_z(const EulerVar& U);

// };



struct VecField {
    
};


using EulerVec = Container1D<double, N_EQS_EULER>; 

struct EulerVecField : public VecField{
    Vector<EulerVec> cell_values;
};


class SolverData{
protected:
    
    unique_ptr<VecField>  solution, solution_old, residual;

public:
    virtual VecField& get_solution() = 0;
    virtual const VecField& get_solution() const = 0;
};



class EulerSolverData : public SolverData{

public:
    EulerSolverData(const Config& config);


    
    
    EulerVecField& get_solution() override {return static_cast<EulerVecField&>(*solution);}

    const EulerVecField& get_solution() const override {return static_cast<EulerVecField&>(*solution);}
};





class EulerEqs{

    static EulerVec prim_to_cons(const EulerVec& V);
    static EulerVec cons_to_prim(const EulerVec& U);

    static double pressure(const EulerVec& U);
    static double sound_speed(const EulerVec& U);

    static EulerVec inviscid_flux_x(const EulerVec& U);
    static EulerVec inviscid_flux_y(const EulerVec& U);
    static EulerVec inviscid_flux_z(const EulerVec& U);
};



