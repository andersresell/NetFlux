#pragma once
#include "../includes.hpp"
#include "StaticContainer.hpp"

/*Dynamic container used to hold solution variables
The type below is used to store vectors of matrices (such as gradients). 
They are stored in this fashion [A0, A1 ... Ai ... AN], where Ai is a MxN matrix. 
M will typically represent the number of solution variables, while  
It is assumed that the number of colums of the matrices is known at compile time*/
template<typename T>
class DynamicContainer3D{
    Index length; //Length of outer vector
    Index M;
    Index N;

    T* data;

public:
    DynamicContainer3D() {};

    DynamicContainer3D(Index length, Index M, Index N) : L{L}, M{M}, N{N} {
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
        return data[l * M * N];
    }

    
    template<typename StaticContainerType>
    void convert_to_static_type(Index i, StaticContainerType& container){
        static_assert(sizeof(T) == 8); //just doubles for now
        assert(sizeof(StaticContainerType) == sizeof(T)*M*N);
        
        container.data

        return static_cast<StaticContainerType*>(data[l * M * N]);
    }

    void operator*=(T rhs){
        for (Index i{0}; i<SIZE; i++) data[i] *= rhs;
    }

    void operator=(T rhs){
        for (Index i{0}; i<SIZE; i++) data[i] = rhs;
    }

    void set_zero() {for (size_t i{0}; i<length*M*N; i++) data[i] = 0;}

    using EigenMat = Eigen::MatrixX<T>;
    /*Maps the matrix data of element l to an eigen matrix*/
    EigenMat& get_eigen_matrix(Index l) {
        Eigen::Map<EigenMat> matrix(data + N*M*l, M, N);
        return matrix;    
    }
    /*Maps the matrix data of element l to an eigen matrix*/
    const EigenMat& get_eigen_matrix(Index l) const {
        Eigen::Map<EigenMat> matrix(data + N*M*l, M, N);
        return matrix;    
    }

    Index length() const{return length;}
    Index rows() const {return M;}
    Index cols() const {return N;}
    

    ~DynamicContainer3D() {delete [] data ;}
};

template<typename T>
class DynamicContainer2D : public DynamicContainer3D<T>{

public:
    DynamicContainer2D() {}


    DynamicContainer2D(Index length, Index M) : DynamicContainer3D(length, M, 1) {}

    T& operator()(Index l, Index i) {
        assert(l<length && i<M);
        return data[l * M + i];
    }

    const T& operator()(Index l, Index i) const {
        assert(l<length && i<M);
        return data[l * M + i];
    }
};