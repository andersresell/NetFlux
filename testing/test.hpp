#pragma once

#include "../include/includes.hpp"

using namespace std;

class Static
{
public:
    double sdata[4] = {0};

    double &operator[](Index i)
    {
        return sdata[i];
    }

    void print()
    {
        for (int i{0}; i < 4; i++)
            cout << sdata[i] << ",";
        cout << std::endl;
    }
};

class Dynamic
{
public:
    double *ddata;

    Dynamic()
    {
        ddata = new double[16];
        for (int i{0}; i < 16; i++)
            ddata[i] = i;
    }

    void print()
    {
        for (int i{0}; i < 16; i++)
            cout << ddata[i] << ",";
        cout << std::endl;
    }

    double* operator[](Index i)
    {
        return &ddata[4 * i];
    }

    //template <typename StaticContainerType>
    // void convert_to_static_type(Index i, Static& s)
    // {
    //     //static_assert(sizeof(StaticContainerType) == 4*sizeof(double));

    //     s = *(Static*)&ddata[i*4];

    // }

    template<typename StaticEigenType>
    Eigen::Map<StaticEigenType> get_variable(Index l){
        // Eigen::Map<StaticEigenType> matrix(ddata + 4*l, 4, 1);
        // return matrix;     
        return Eigen::Map<StaticEigenType>(ddata + 4*l);
        //return Eigen::Map<StaticEigenType>(ddata + 4*l, 4, 1);
    }
};