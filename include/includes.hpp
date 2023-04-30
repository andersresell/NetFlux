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



using Index = uint32_t; //Used for accessing cell indices, etc 
using ShortIndex = uint16_t; //Used for looping over shorter numbers as spatial dimensions, etc

using string = std::string;
using Vec3 = Eigen::Vector3d;
inline string horizontal_string_Vec3(const Vec3& v) {return "["+std::to_string(v[0])+", "+std::to_string(v[1])+", "+std::to_string(v[2])+" ]";}
template <typename T, size_t N> using array = std::array<T,N>;
template <typename TA, typename TB> using map = std::map<TA,TB>;
using std::cout;
using std::endl;

#define FAIL_IF(EXP) ({if (EXP) exit(EXIT_FAILURE);})
#define FAIL_IF_MSG(EXP, MSG) ({if (EXP) {std::cerr << MSG << "\n"; exit(EXIT_FAILURE);}})
#define FAIL_MSG(MSG) ({std::cerr << MSG << "\n"; exit(EXIT_FAILURE);})



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
    assert (a.size() == b.size());
    for (size_t i{0}; i<a.size(); i++) 
        if (a[i] != b[i]) return false;
    return true;
}
