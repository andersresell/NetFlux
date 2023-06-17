#pragma once

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <map>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <set>
#include <utility>
#include <stdint.h>
#include <memory>
#include <any>
#include <chrono>
#include <cfloat>

using std::unique_ptr;
using std::make_unique;
using std::string;
using std::array;
using std::map;
using std::move;
using std::cout;
using std::endl;
using std::max;
using std::min;


using Index = uint32_t; //Used for accessing cell indices, etc 
using ShortIndex = uint16_t; //Used for looping over shorter numbers as spatial dimensions, etc

using Vec3 = Eigen::Vector3d;

inline string horizontal_string_Vec3(const Vec3& v) {return "["+std::to_string(v[0])+", "+std::to_string(v[1])+", "+std::to_string(v[2])+" ]";}

constexpr ShortIndex N_DIM{3}; //spatial dimensions

#define FAIL_IF(EXP) ({if (EXP) exit(EXIT_FAILURE);})
#define FAIL_IF_MSG(EXP, MSG) ({if (EXP) {std::cerr << MSG << "\n"; exit(EXIT_FAILURE);}})
#define FAIL_MSG(MSG) ({std::cerr << MSG << "\n"; exit(EXIT_FAILURE);})
#define assert_msg(EXP, MSG) ({if (!EXP) {std::cerr << MSG << "\n"; assert(false);}})


/*Adding range checking to the [] operator for std::vector in debug mode (NDEBUG not defined)*/
template <typename T>
class Vector final : public std::vector<T> {
public:
    using std::vector<T>::vector;
    
    #ifndef NDEBUG

    T& operator[](size_t i) {
        return this->at(i);
    }
    const T& operator[](size_t i) const {
        return this->at(i);
    }

    #endif
};


template<typename T, size_t N>
inline bool arrays_equal(const array<T, N>& a,const array<T, N>& b){
    for (size_t i{0}; i<N; i++) 
        if (a[i] != b[i]) return false;
    return true;
}

template<typename T>
inline int sign(T val){
    return val >= 0;
}