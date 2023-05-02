#pragma once

#include <iostream>
#include <vector>
#include <assert.h>
#include <memory>
using Index = size_t;
constexpr Index N_TET_NODES=4;

template<typename T>
using Vector = std::vector<T>;

using namespace std;


template<typename T, Index N>
class Container1D{};

template<typename T, Index M, Index N>
class Container2D{};

template<Index N>
class FlowVec : public Container1D<N> {};

template<Index N>
class FlowGrad : public: Container2D<3,N> {};

template<Index N>
class FlowVecField{
    Vector<FlowVec<N>> vec_field;
};

template<Index N>
class FlowGradField{
    Vector<FlowGrad<N>> grad_field;
};

class EulerVecField : public FlowVecField<5>{
    using EulerVec = FlowVec<5>;
public:
    static void get_inv_flux(const EulerVec& U);
};

class NS_VecField : public EulerVecField{
    using NS_Vec = FlowVec<5>;
    void do_something(const NS_Vec& U){ get_inv_flux(U);}
};

class SolverBase{

};

template<Index N>
class Solver : SolverBase{
protected:
    unique_ptr<FlowVecField<N>> solution, solution_old;
    virtual void func_something() = 0;
};

class EulerSolver : Solver<5>{
    EulerSolver(){
        solution = make_unique<EulerVecField>();
    }
};


class Driver{
    unique_ptr<SolverBase> solver;

};






// opt 2:
template<typename N>
using FlowVar = Container1D<double, N>;
using EulerVar = FlowVar<5>;
using NS_Var = EulerVar;
struct EulerEquations{

    static EulerVar convective_flux(EulerVar U);
};

struct NS_Equations : public EulerEquations{
    static double get_viscosity(NS_Var U); 
};

struct BaseSolver{

};

template <typename T>
struct Solver2 : public BaseSolver{
    Solver2() {}
    unique_ptr<T> solution, solution_old;
};


struct EulerSolver : public Solver2<EulerVar>{
    EulerSolver() : Solver2() {solution = make_unique<Vector<EulerVar>>();}

};

struct NS_Solver : public EulerSolver{

};