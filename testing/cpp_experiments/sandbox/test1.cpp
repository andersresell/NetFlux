#include <iostream>
#include <vector>
// #include <sstream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <eigen3/Eigen/Dense>
#include <mpi.h>
using namespace std;
// Define a serialization function for Eigen::Vector3d

using Scalar = double;

using namespace std;

template <typename T>
static MPI_Datatype get_MPI_Datatype()
{
    using std::is_same;
    if constexpr (is_same<T, double>::value)
        return MPI_DOUBLE;
    else if constexpr (is_same<T, uint32_t>::value)
        return MPI_UINT32_T;
    else if constexpr (is_same<T, char>::value)
        return MPI_BYTE;
    else
    {
        static_assert(false);
    }
}

int main()
{
    cout << "hello\n";
}