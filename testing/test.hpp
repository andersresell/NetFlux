#pragma once

#include <iostream>
#include <vector>
#include <assert.h>

using Index = size_t;
constexpr Index N_TET_NODES=4;

template<typename T>
using Vector = std::vector<T>;

using namespace std;


struct Face{
    Index i,j; //Indices of cell i and j
};

struct Cell{
    double cell_volume;
};
