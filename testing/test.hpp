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


template <typename T>
class BaseClass {
public:
    template <typename U>
    virtual T getMemberX() = 0;
};

class DerivedClass : public BaseClass<int> {
public:
    template <typename U>
    int getMemberX() override {
        // implementation
    }
};
