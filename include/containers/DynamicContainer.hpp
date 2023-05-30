#pragma once
#include "../includes.hpp"
#include "StaticContainer.hpp"

/*Dynamic container used to hold solution variables
The type below is used to store vectors of matrices (such as gradients). 
They are stored in this fashion [A0, A1 ... Ai ... AN], where Ai is a MxN matrix. 
M will typically represent the number of solution variables, while  
It is assumed that the number of colums of the matrices N is known at compile time*/

template<typename T, ShortIndex N>
class DynamicContainer3D{

    Index length; //Length of outer vector
    Index M; //number of rows of each variable, 

    T* data;

public:
    DynamicContainer3D() {};

    DynamicContainer3D(Index length, Index M) : L{L}, M{M} {
        data = new T[length * M * N](0);
    }

    T& operator()(Index l, Index i, Index j) {
        assert(l<length && i<M && j<N);
        return data[l * M*N + i*N + j];
    }

    const T& operator()(Index l, Index i, Index j) const {
        assert(l<length && i<M && j<N);
        return data[l * M*N + i*N + j];
    }

    /*Returns pointer to matrix l*/
    T* operator[](Index l) {
        assert(l < length);
        return data + l * M * N;
    }
    const T* operator[](Index l) const {
        assert(l < length);
        return data + l * M * N;
    }

    void operator*=(T rhs){
        for (Index i{0}; i<SIZE; i++) data[i] *= rhs;
    }

    void operator=(T rhs){
        for (Index i{0}; i<SIZE; i++) data[i] = rhs;
    }

    void set_zero() {for (size_t i{0}; i<length*M*N; i++) data[i] = 0;}


    template<typename StaticEigenType>
    void set_variable(Index l, const StaticEigenType& variable){
        assert(StaticEigenType::RowsAtCompileTime == M && StaticEigenType::ColsAtCompileTime == N);
        for (Index i{0}; i<M; i++)
            for (Index j{0}; j<N; j++)
                data[l * M*N + i*N + j] = variable(i,j);

    }

    template<typename StaticEigenType>
    Eigen::Map<StaticEigenType> get_variable(Index l){
        assert(StaticEigenType::RowsAtCompileTime == M && StaticEigenType::ColsAtCompileTime == N);
        return Eigen::Map<StaticEigenType>(data + N*M*l);
    }

    template<typename StaticEigenType>
    const Eigen::Map<StaticEigenType> get_variable(Index l) const{
        assert(StaticEigenType::RowsAtCompileTime == M && StaticEigenType::ColsAtCompileTime == N);
        return Eigen::Map<StaticEigenType>(data + N*M*l);
    }




    Index length() const{return length;}
    Index rows() const {return M;}
    Index cols() const {return N;}
    

    ~DynamicContainer3D() {delete [] data ;}
};

template<typename T>
class DynamicContainer2D : public DynamicContainer3D<T, 1>{

public:
    DynamicContainer2D() {}


    DynamicContainer2D(Index length, Index M) : DynamicContainer3D(length, M) {}

    T& operator()(Index l, Index i) {
        assert(l<length && i<M);
        return data[l * M + i];
    }

    const T& operator()(Index l, Index i) const {
        assert(l<length && i<M);
        return data[l * M + i];
    }
};