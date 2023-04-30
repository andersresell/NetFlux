#pragma once

#include <iostream>
#include <vector>
#include <assert.h>

using Index = size_t;
constexpr Index N_TET_NODES=4;

template<typename T>
using Vector = std::vector<T>;

template<typename T, Index N>
class Container{

protected:
    T data[N];
public:
    Container() : data{} {}

    Container(std::initializer_list<T> init){
        assert(init.size() == N);
        std::copy(init.begin(), init.end(), data);
    }
    
    T& operator[](Index i) {
        assert(i<N);
        return data[i];
    }
    const T& operator[](Index i) const{
        assert(i<N);
        return data[i];
    }
    
    Container<T,N> operator+(Container<T,N> rhs) const {
        for (Index i{0}; i<N; i++) rhs.data[i] += data[i];
        return rhs;
    };

    Container<T,N> operator-(const Container<T,N>& rhs) const {
        Container lhs{*this};
        for (Index i{0}; i<N; i++) lhs.data[i] -= rhs.data[i];
        return lhs;
    };

    void operator-=(const Container<T,N>& rhs) {
        for (Index i{0}; i<N; i++) data[i] -= rhs.data[i];
    };

    void operator+=(const Container<T,N>& rhs) {
        for (Index i{0}; i<N; i++) data[i] += rhs.data[i];
    };

    Container<T,N> operator*(T rhs) const {
        Container lhs{*this};
        for (Index i{0}; i<N; i++) lhs.data[i] *= rhs;
        return lhs;
    };

    void operator*=(T rhs){
        for (Index i{0}; i<N; i++) data[i] *= rhs;
    }

    Index size() const {return N;}
    
    friend std::ostream& operator<<(std::ostream& os, const Container& rhs){
        for (Index i{0}; i<N; i++) os << rhs[i] << " ";
        return os << std::endl;
    }
};




template<Index N_EQS>
class Equation{
public:
    struct Variable : public Container<double, N_EQS>{

    };    
};



constexpr Index N_EQS_EULER{5};
class EulerEquation : public Equation<N_EQS_EULER>{
public:
    struct EulerVariable : public Variable{
        double get_density() const {return data[0];}
    };
};

class NS_EQUATION : public EulerEquation{
    struct NS_VARIABLE : public EulerVariable{
        double get_visc_or_something(){return 0.3123113;}
    };
};


template<typename Equation>
class Solver{
public:
    Vector<typename Equation::Variable> solution;


};

class EulerSolver : public Solver<typename EulerEquation>{
    double get_density(Index i) {return this->solution[i].g}
};