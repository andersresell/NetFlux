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
using VecType = Container1D<double, N_EQS>;

template<ShortIndex N_EQS>
using GradType = Container2D<double, N_DIM, N_EQS>;


struct VecField {};

struct GradField {};


using EulerVec = VecType<N_EQS_EULER>; 

struct EulerVecField : public VecField{
    Vector<EulerVec> cell_values;
    EulerVecField();
    //Vector<EulerVec>* get_cell_values() {return &cell_values;} 
};

using EulerGrad = GradType<N_EQS_EULER>;

struct EulerGradField : public GradField{
    Vector<EulerGrad> cell_values;
};


class SolverData{
protected:

    unique_ptr<VecField>  solution, solution_old, residual, primvars;
    unique_ptr<GradField> primvars_gradient;

    SolverData() = default;
public:

    virtual VecField& get_solution() = 0;
    virtual const VecField& get_solution() const = 0;
};



class EulerSolverData : public SolverData{

public:
    EulerSolverData(const Config& config);


    
    
    EulerVecField& get_solution() override {return static_cast<EulerVecField&>(*solution);}

    const EulerVecField& get_solution() const override {
        return static_cast<EulerVecField&>(*solution);}
};





namespace EulerEqs{

    EulerVec prim_to_cons(const EulerVec& V);
    EulerVec cons_to_prim(const EulerVec& U);

    double pressure(const EulerVec& U);
    double sound_speed(const EulerVec& U);

    double projected_velocity(const EulerVec& U, const Vec3& normal);


    /*returns |V| + c where V = <velocity, normal> and c = sound speed */
    double conv_spec_rad(const EulerVec& U, const Vec3& normal);

    EulerVec inviscid_flux(const EulerVec& U, const Vec3& normal);


    
    inline EulerVec prim_to_cons(const EulerVec& V){
        return {V[0], 
                V[0]*V[1],
                V[0]*V[2],
                V[0]*V[3],
                GAMMA_MINUS_ONE_INV * V[4] + 0.5*V[0]*(V[1]*V[1] + V[2]*V[2] + V[3]*V[3])};
    }
    inline EulerVec cons_to_prim(const EulerVec& U){
        return {U[0], 
                U[1]/U[0],
                U[2]/U[0],
                U[3]/U[0],
                pressure(U)};
    }    


    inline double pressure(const EulerVec& U){
        return GAMMA_MINUS_ONE * (U[4] - 0.5/U[0]*(U[1]*U[1] + U[2]*U[2] + U[3]*U[3]));
    }

    inline double sound_speed(const EulerVec& U){
        return sqrt(GAMMA * pressure(U) / U[0]);
    }

    inline double projected_velocity(const EulerVec& U, const Vec3& normal){
        return (U[1]*normal.x() + U[2]*normal.y() + U[3]*normal.z())/U[0];
    }

    inline double conv_spectral_rad(const EulerVec& U, const Vec3& normal){
        return abs(projected_velocity(U, normal)) + sound_speed(U);
    }

    inline EulerVec inviscid_flux(const EulerVec& U, const Vec3& normal){
        double rho = U[0];
        double p = pressure(U);
        double V_normal = (U[1]*normal.x() + U[2]*normal.y() + U[3]*normal.z())/rho;

        return {V_normal*rho,
                V_normal*U[1] + p*normal.x(),
                V_normal*U[2] + p*normal.y(),
                V_normal*U[3] + p*normal.z(),
                V_normal*(U[4] + p)};
    }


}
