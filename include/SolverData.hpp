#pragma once
#include "includes.hpp"
#include "containers/DynamicContainer.hpp"
#include "containers/StaticContainer.hpp"
#include "Config.hpp"






struct VecField final : public DynamicContainer2D<double>{
    VecField(Index N_CELLS, ShortIndex N_EQS) : DynamicContainer2D(N_CELLS, N_EQS) {}
    ShortIndex get_N_EQS() const {return rows();}
};

struct GradField final : public DynamicContainer3D<double, N_DIM>{
    GradField(Index N_CELLS, ShortIndex N_EQS) : DynamicContainer3D(N_CELLS, N_EQS) {}
    ShortIndex get_N_EQS() const {return rows();}
}; 




class SolverData{
protected:

    unique_ptr<VecField>  solution, 
                          solution_old, 
                          flux_balance, 
                          primvars;

    unique_ptr<GradField> primvars_gradient;

    unique_ptr<VecField> primvars_limiter,
                         primvars_max, 
                         primvars_min;


    SolverData() = default;
public:

    VecField& get_solution() {return *solution;}
    const VecField& get_solution() const {return *solution;}

    VecField& get_flux_balance() {return *flux_balance;}
    const VecField& get_flux_balance() const {return *flux_balance;}

    VecField& get_primvars() {return *primvars;}
    const VecField& get_primvars() const {return *primvars;}

    GradField& get_primvars_gradient() {return *primvars_gradient;}
    const GradField& get_primvars_gradient() const {return *primvars_gradient;}

    VecField& get_primvars_limiter() {return *primvars_limiter;}
    const VecField& get_primvars_limiter() const {return *primvars_limiter;}

    VecField& get_primvars_max() {return *primvars_max;}
    const VecField& get_primvars_max() const {return *primvars_max;}

    VecField& get_primvars_min() {return *primvars_min;}
    const VecField& get_primvars_min() const {return *primvars_min;}

    virtual SolverType get_solver_type() const = 0; 
};



template<ShortIndex N_EQS>
using FlowVec = Eigen::Vector<double, N_EQS>;
template<ShortIndex N_EQS>
using FlowGrad = Eigen::Matrix<double, N_EQS, N_DIM>;

// using EulerVec = StaticContainer1D<double, N_EQS_EULER>; 
// using EulerGrad = StaticContainer2D<double, N_EQS_EULER, N_DIM>; 
using EulerVec = FlowVec<N_EQS_EULER>;
using EulerGrad = FlowGrad<N_EQS_EULER>;

using EulerVecMap = Eigen::Map<EulerVec>;
using EulerGradMap = Eigen::Map<EulerGrad>;

class EulerSolverData : public SolverData{

    EulerVec U_L, U_R, V_L, V_R;
    EulerVec Flux_inv;

    Vector<Vec3> Delta_S; //Used in time step computation
    

public:
    EulerSolverData(const Config& config);

    virtual SolverType get_solver_type() const {return SolverType::Euler;}

    EulerVecMap get_U_L_map() {return EulerVecMap(U_L.data());}
    EulerVecMap get_U_R_map() {return EulerVecMap(U_R.data());}
    EulerVecMap get_V_L_map() {return EulerVecMap(V_L.data());}
    EulerVecMap get_V_R_map() {return EulerVecMap(V_R.data());}
    EulerVecMap get_Flux_inv_map() {return EulerVecMap(Flux_inv.data());}

    EulerVecMap get_empty_map() const {return EulerVecMap(nullptr);}

    Vector<Vec3>& get_Delta_S() {return Delta_S;}
    const Vector<Vec3>& get_Delta_S() const {return Delta_S;}
};

    /*Discontinuing StaticContainer, using Eigen instead*/

    // template<ShortIndex N_EQS>
    // using FlowVec = StaticContainer1D<double, N_EQS>;
    // template<ShortIndex N_EQS>
    // using FlowGrad = StaticContainer2D<double, N_EQS, N_DIM>;




namespace EulerEqs{

    constexpr double GAMMA{1.4};
    constexpr double GAMMA_MINUS_ONE{1-GAMMA};
    constexpr double GAMMA_MINUS_ONE_INV{1/GAMMA_MINUS_ONE};

    /*Templating the various functions since the EulerVecType can be either EulerVec or EulerVecMap*/

    template<typename EulerVecType>
    inline void prim_to_cons(const EulerVecType& V, EulerVecType& U);

    template<typename EulerVecType>
    inline void cons_to_prim(const EulerVec& U, EulerVecType& V);


    template<typename EulerVecType>
    inline double pressure(const EulerVecType& U);

    template<typename EulerVecType>
    inline double sound_speed_conservative(const EulerVecType& U);


    template<typename EulerVecType>
    inline double sound_speed_primitive(const EulerVecType& V);


    template<typename EulerVecType>
    inline double projected_velocity(const EulerVecType& U, const Vec3& normal);


    /*returns |V| + c where V = <velocity, normal> and c = sound speed */
    template<typename EulerVecType>
    inline double conv_spec_rad(const EulerVecType& U, const Vec3& normal);

    template<typename EulerVecType>
    inline EulerVec inviscid_flux(const EulerVecType& U, const Vec3& normal);


    template<typename EulerVecType>
    inline void prim_to_cons(const EulerVecType& V, EulerVecType& U){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER and EulerVecType::ColsAtCompileTime == 1);
        U = {V[0], 
            V[0]*V[1],
            V[0]*V[2],
            V[0]*V[3],
            GAMMA_MINUS_ONE_INV * V[4] + 0.5*V[0]*(V[1]*V[1] + V[2]*V[2] + V[3]*V[3])};
        // U[0] = V[0];
        // U[1] = V[0]*V[1];
        // U[2] = V[0]*V[2],
        // U[3] = V[0]*V[3],
        // U[4] = GAMMA_MINUS_ONE_INV * V[4] + 0.5*V[0]*(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
    }

    template<typename EulerVecType>
    inline void cons_to_prim(const EulerVecType& U, EulerVecType& V){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER and EulerVecType::ColsAtCompileTime == 1);
        V = {U[0], 
                U[1]/U[0],
                U[2]/U[0],
                U[3]/U[0],
                pressure(U)};
        // V[0] = U[0], 
        // V[1] = U[1]/U[0],
        // V[2] = U[2]/U[0],
        // V[3] = U[3]/U[0],
        // V[4] = pressure(U)};
    }    


    template<typename EulerVecType>
    inline double pressure(const EulerVecType& U){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        return GAMMA_MINUS_ONE * (U[4] - 0.5/U[0]*(U[1]*U[1] + U[2]*U[2] + U[3]*U[3]));
    }


    template<typename EulerVecType>
    inline double sound_speed_conservative(const EulerVecType& U){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        return sqrt(GAMMA * pressure(U) / U[0]);
    }

    template<typename EulerVecType>
    inline double sound_speed_primitive(const EulerVecType& V){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        return sqrt(GAMMA * V[4] / V[0]);
    }

    template<typename EulerVecType>
    inline double projected_velocity(const EulerVecType& U, const Vec3& normal){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        return (U[1]*normal.x() + U[2]*normal.y() + U[3]*normal.z())/U[0];
    }


    template<typename EulerVecType>
    inline double conv_spectral_rad(const EulerVecType& U, const Vec3& normal){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);
        return abs(projected_velocity(U, normal)) + sound_speed(U);
    }

    template<typename EulerVecType>
    inline EulerVec inviscid_flux(const EulerVecType& U, const Vec3& normal, EulerVecType& Flux){
        static_assert(EulerVecType::RowsAtCompileTime == N_EQS_EULER && EulerVecType::ColsAtCompileTime == 1);

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
