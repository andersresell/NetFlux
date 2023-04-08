
#include "../include/Utilities.hpp"
#include <iostream>

int main(){
    FlowVec5 f{1,2,3,4,5};
    std::cout << "hello fucker\n";
    std::cout << f.u1<<std::endl;

    Eigen::Vector3d v;
    v << 1,2,3;
    std::cout << v <<std::endl;
}