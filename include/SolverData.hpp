#pragma once
#include "includes.hpp"
#include "containers/DynamicContainer.hpp"
#include "containers/StaticContainer.hpp"
#include "Config.hpp"






struct VecField final : public DynamicContainer2D<double>{
    VecField(Index N_CELLS, ShortIndex N_EQS) : DynamicContainer2D(N_CELLS, N_EQS) {}
    ShortIndex get_N_EQS() const {return rows();}
};

struct GradField final : public DynamicContainer3D<double>{
    GradField(Index N_CELLS, ShortIndex N_EQS) : DynamicContainer3D(N_CELLS, N_EQS, N_DIM) {}
    ShortIndex get_N_EQS() const {return rows();}
}; 




class SolverData{
protected:

    unique_ptr<VecField>  solution, solution_old, residual, primvars;
    unique_ptr<GradField> primvars_gradient;

    SolverData() = default;
public:

    VecField& get_solution() {return *solution;}
    const VecField& get_solution() const {return *solution;}

    VecField& get_residual() {return *residual;}
    const VecField& get_residual() const {return *residual;}

    VecField& get_primvars() {return *primvars;}
    const VecField& get_primvars() const {return *primvars;}

    GradField& get_primvars_gradient() {return *primvars_gradient;}
    const GradField& get_primvars_gradient() const {return *primvars_gradient;}

};



class EulerSolverData : public SolverData{

public:
    EulerSolverData(const Config& config);

};




    using EulerVec = StaticStaticContainer1D<double, N_EQS_EULER>; 
    using EulerGrad = StaticContainer2D<double, N_EQS_EULER, N_DIM>; 

namespace EulerEqs{

    constexpr double GAMMA{1.4};
    constexpr double GAMMA_MINUS_ONE{1-GAMMA};
    constexpr double GAMMA_MINUS_ONE_INV{1/GAMMA_MINUS_ONE};


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
