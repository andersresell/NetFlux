#pragma once

#include <iostream>
#include <vector>
#include <assert.h>

using Index = size_t;
constexpr Index N_TET_NODES=4;

template <typename T>
class Vector final : public std::vector<T> {
public:
    using std::vector<T>::vector;

    #ifndef NDEBUG

    T& operator[](size_t i) {
        std::cout << "DEBUG VERSION\n";
        return this->at(i);
    }

    #endif
    


};

template<typename T, size_t N>
struct Container{

    T data[N];
    Container() : data{} {}

    Container(std::initializer_list<T> init){
        assert(init.size() == N);
        std::copy(init.begin(), init.end(), data);
    }

    
    T& operator[](size_t i) {
        assert(i<N);
        return data[i];
    }
    const T& operator[](size_t i) const{
        assert(i<N);
        return data[i];
    }

    
    Container<T,N> operator+(Container<T,N> rhs) const {
        for (size_t i{0}; i<N; i++) rhs.data[i] += data[i];
        return rhs;
    };

    Container<T,N> operator-(const Container<T,N>& rhs) const {
        Container lhs{*this};
        for (size_t i{0}; i<N; i++) lhs.data[i] -= rhs.data[i];
        return lhs;
    };

    void operator-=(const Container<T,N>& rhs) {
        for (size_t i{0}; i<N; i++) data[i] -= rhs.data[i];
    };

    void operator+=(const Container<T,N>& rhs) {
        for (size_t i{0}; i<N; i++) data[i] += rhs.data[i];
    };

    Container<T,N> operator*(T rhs) const {
        Container lhs{*this};
        for (size_t i{0}; i<N; i++) lhs.data[i] *= rhs;
        return lhs;
    };

    void operator*=(T rhs){
        for (size_t i{0}; i<N; i++) data[i] *= rhs;
    }

    size_t size() const {return N;}
    
    friend std::ostream& operator<<(std::ostream& os, const Container& rhs){
        for (size_t i{0}; i<N; i++) os << rhs[i] << " ";
        return os << std::endl;
    }
};


struct TetConnect : public Container<Index, N_TET_NODES> {
    
    TetConnect(std::initializer_list<Index> init) : Container{init} {}
    Index& a() {return data[0];}
    Index& b() {return data[1];}
    Index& c() {return data[2];}
    Index& d() {return data[3];}
};