#include <iostream>
#include <vector>
// #include <sstream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <eigen3/Eigen/Dense>
#include <array>
using namespace std;
// Define a serialization function for Eigen::Vector3d

using Scalar = double;

using namespace std;

int main()
{
    cout << "hello\n";

    // vector<int> v{1, 2, 3};
    // v[5] = 4;

    array<int, 3> a{1, 2, 3};
    a[4] = 5;
}