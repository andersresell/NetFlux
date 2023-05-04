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
template<ShortIndex N_EQS>
using FlowVar = Container1D<double, N_EQS>;

template<ShortIndex N_EQS>
using FlowGrad = Container2D<double, N_DIM, N_EQS>;





template<ShortIndex N_EQS> 
struct FlowField{
    Vector<FlowVar<N_EQS>> cell_values;    


    template <typename T>
    virtual T& getMemberX() = 0;
};

template<ShortIndex N_EQS>
struct FlowGradField{
    Vector<FlowGrad<N_EQS>> cell_val_gradiens;
};


using EulerVar = FlowVar<N_EQS_EULER>;

struct EulerField : public FlowField<N_EQS_EULER>{

    EulerField(Config & config);

    static EulerVar prim_to_cons(const EulerVar& V);
    static EulerVar cons_to_prim(const EulerVar& U);

    static double pressure(const EulerVar& U);
    static double sound_speed(const EulerVar& U);

    static EulerVar inviscid_flux_x(const EulerVar& U);
    static EulerVar inviscid_flux_y(const EulerVar& U);
    static EulerVar inviscid_flux_z(const EulerVar& U);
};



/*perhaps for later*/
struct NS_Field : public EulerField{
    
};
