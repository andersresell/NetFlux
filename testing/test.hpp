#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <assert.h>
#include <memory>
using Index = size_t;
constexpr Index N_TET_NODES=4;

template<typename T>
using Vector = std::vector<T>;

using namespace std;

class A{
public:
    int a,b;
    double d;
    string str;
    A() = default;

};

class MyClass {
public:
    int value;
    A a;

    // User-defined constructor
    MyClass(int val) {
        value = val;
        std::cout << "User-defined constructor called" << std::endl;
    }

    // Default constructor explicitly requested
    MyClass() = default;
};