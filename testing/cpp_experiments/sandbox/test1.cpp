#include <iostream>
#include <vector>
// #include <sstream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <eigen3/Eigen/Dense>
using namespace std;
// Define a serialization function for Eigen::Vector3d

using Scalar = double;

using namespace std;

template <size_t N>
class B
{
public:
    size_t data{3};
    size_t get_data() { return data; }
};

class D : public B<3>
{

    void printdata() { cout << "data: " << data; }
};

template <size_t N>
void func(B<N> baseclass)
{
    cout << "func " << baseclass.get_data() << endl;
}

int main()
{
    D d;
    func(d);
}