#pragma once

#include "../include/includes.hpp"

using namespace std;

template<typename T, size_t N>
class A{

    T* data;

public:

    A(){data = new T[N];}

    void operator=(T rhs){
        for (int i{0}; i<N; i++) data[i] = rhs;
    }
    void operator+=(T rhs){
                for (int i{0}; i<N; i++) data[i] += rhs;
    }

    void print() {for (int i{0}; i<N; i++) cout << data[i] << ", "; cout << endl;}
};


class B : public A<double, 10>{

public:
using A<double, 10>::operator=; // Add this line to bring the base class operator= into scope
    B() : A<double,10>() {}
};