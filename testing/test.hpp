#pragma once

#include <iostream>
#include <vector>
#include <assert.h>

using Index = size_t;
constexpr Index N_TET_NODES=4;

template<typename T>
using Vector = std::vector<T>;

using namespace std;

template<typename T, Index N>
class Container{



protected:
    T data[N];
public:
    Container() : data{} {cout << "constructing Container of zeros\n";}


    template<typename... Args>
    Container(Args... args) : data{static_cast<T>(args)...} {
        static_assert(sizeof...(args) == N);
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


constexpr Index N_EQS_EULER{5};

struct EulerVar : public Container<double, N_EQS_EULER>{
    EulerVar() =  default;
    //EulerVar(std::initializer_list<double> init) : Container(init) {cout << "constructing EulerVar init list\n";}
    
    template<typename... Args>
    
    EulerVar(Args&&... args) : Container(std::forward<Args>(args)...) {}

    static EulerVar prim_to_cons(const EulerVar& V);
    static EulerVar cons_to_prim(const EulerVar& U);

    static double pressure(const EulerVar& U);
    static double sound_speed(const EulerVar& U);

    static EulerVar inviscid_flux_x(const EulerVar& U);
    static EulerVar inviscid_flux_y(const EulerVar& U);
    static EulerVar inviscid_flux_z(const EulerVar& U);

};