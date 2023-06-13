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




template<typename T, size_t cols_>
class DynamicContainer3D{

    size_t size_; //Length of outer vector
    size_t rows_; //number of rows of each variable, 

    T* data;

public:

    DynamicContainer3D(size_t size, size_t rows) : size_{size}, rows_{rows} {
        data = new T[size_ * rows_ * cols_](0);
    }


    template<typename StaticEigenType>
    Eigen::Map<StaticEigenType> get_variable(size_t i){
        assert(StaticEigenType::RowsAtCompileTime == rows_ && StaticEigenType::ColsAtCompileTime == cols_);
        return Eigen::Map<StaticEigenType>(data + rows_ * cols_ * i);
    }



    ~DynamicContainer3D() {delete [] data ;}
};

template<typename T>
class DynamicContainer2D : public DynamicContainer3D<T, 1>{
public:
    DynamicContainer2D(size_t size, size_t rows) : DynamicContainer3D<T,1>(size, rows) {};
};


// template <Index size>
// class View{
//     double* data;

// public:
//     View(double* data) {data = data;}

//     void operator*=(double rhs) {
//         for (int i{0};i<size;i++) data[i] *= rhs;
//     }

//     View(const View& rhs){
//         data = double[size];
//         std::copy(rhs.data, rhs.data + size, data);
//     }


// };

// class Container{
//     double* data;

//     template <Index size>
//     View<size> get_variable(Index i){
//         return View{data + i*size};
//     }
// };