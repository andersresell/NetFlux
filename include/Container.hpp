#pragma once

#include "includes.hpp"


/*Generic vector type used directly or for building other data types*/
template<typename T, Index N>
class Container1D{

protected:
    T data[N];
public:

    Container1D() : data{} {}

    template<typename... Args>
    Container1D(Args... args) : data{static_cast<T>(args)...} {
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
    
    Container1D<T,N> operator+(Container1D<T,N> rhs) const {
        for (Index i{0}; i<N; i++) rhs.data[i] += data[i];
        return rhs;
    };

    Container1D<T,N> operator-(const Container1D<T,N>& rhs) const {
        Container1D lhs{*this};
        for (Index i{0}; i<N; i++) lhs.data[i] -= rhs.data[i];
        return lhs;
    };

    void operator-=(const Container1D<T,N>& rhs) {
        for (Index i{0}; i<N; i++) data[i] -= rhs.data[i];
    };

    void operator+=(const Container1D<T,N>& rhs) {
        for (Index i{0}; i<N; i++) data[i] += rhs.data[i];
    };

    Container1D<T,N> operator*(T rhs) const {
        Container1D lhs{*this};
        for (Index i{0}; i<N; i++) lhs.data[i] *= rhs;
        return lhs;
    };

    void operator*=(T rhs){
        for (Index i{0}; i<N; i++) data[i] *= rhs;
    }

    Index size() const {return N;}
    
    friend std::ostream& operator<<(std::ostream& os, const Container1D& rhs){
        for (Index i{0}; i<N; i++) os << rhs[i] << " ";
        return os << std::endl;
    }
};


/*Generic matrix type used directly or for building other data types*/
template<typename T, Index M, Index N>
class Container2D{
    constexpr size_t SIZE = M * N;
protected:
    T data[M * N];
public:
    Container2D() : data{} {}

    template<typename... Args>
    Container2D(Args... args) : data{args...} {}

    template<typename... Args>
    Container2D(Args... args) : data{static_cast<T>(args)...} {
        static_assert(sizeof...(args) == SIZE);
    }
    
    T& operator()(Index i, Index j) {
        assert(i<M && j<N);
        return data[i*N+j];
    }
    const T& operator()(Index i, Index j) {
        assert(i<M && j<N);
        return data[i*N+j];
    }
    
    Container1D<T,N> operator+(Container1D<T,N> rhs) const {
        for (Index i{0}; i<SIZE; i++) rhs.data[i] += data[i];
        return rhs;
    };

    Container1D<T,N> operator-(const Container1D<T,N>& rhs) const {
        Container1D lhs{*this};
        for (Index i{0}; i<SIZE; i++) lhs.data[i] -= rhs.data[i];
        return lhs;
    };

    void operator-=(const Container1D<T,N>& rhs) {
        for (Index i{0}; i<SIZE; i++) data[i] -= rhs.data[i];
    };

    void operator+=(const Container1D<T,N>& rhs) {
        for (Index i{0}; i<SIZE; i++) data[i] += rhs.data[i];
    };

    // Container1D<T,N> operator*(T rhs) const {
    //     Container1D lhs{*this};
    //     for (Index i{0}; i<N; i++) lhs.data[i] *= rhs;
    //     return lhs;
    // };

    void operator*=(T rhs){
        for (Index i{0}; i<SIZE; i++) data[i] *= rhs;
    }

    Index size() const {return M*N;}
    
    friend std::ostream& operator<<(std::ostream& os, const Container2D& rhs){
        for (Index i{0}; i<N; i++){ 
            for (Index j{0}; i<M; j++) 
                os << rhs(i,j) << " ";
            os << "\n";
        }
        return os << "\n";
    }
};
