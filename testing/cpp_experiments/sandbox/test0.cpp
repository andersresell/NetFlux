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
enum class HEHE
{
    ja,
    nei,
    joa
};

class PrimalGrid
{
    friend class boost::serialization::access;

public:
    std::vector<Eigen::Vector3d> nodes;

public:
    PrimalGrid()
    {

        nodes.push_back(Eigen::Vector3d(1.0, 2.0, 3.0));
        nodes.push_back(Eigen::Vector3d(4.0, 5.0, 6.0));
        nodes.push_back(Eigen::Vector3d(7.0, 8.0, 9.0));
    }
    void print()
    {
        cout << "printing" << endl;
        for (auto &n : nodes)
            cout << n << endl;
    }
    // Define the serialization function for PrimalGrid

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) { ar &nodes; };
};

namespace boost
{
    namespace serialization
    {

        template <class Archive>
        void serialize(Archive &ar, Eigen::Vector3d &v, const unsigned int version)
        {
            ar &v[0];
            ar &v[1];
            ar &v[2];
        }

    } // namespace serialization
} // namespace boost

using namespace std;
int main()
{
    Array<int, 3> arr{1, 2, 3};
    int a = 1;
    for (size_t i{0}; i < 10; i++)
        arr[i] = i * i;

    for (size_t i{0}; i < 5; i++)
        cout << arr[i] << endl;
    cout << a << endl;
    for (int e : arr)
        cout << e << endl;

    cout << "halla\n";
    // Create a PrimalGrid and add some Eigen::Vector3d objects
    PrimalGrid grid{};
    cout << "first_print\n";
    grid.print();
    // Serialize PrimalGrid to bytes
    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << grid;
    grid.nodes[0] = {3000, 3000, 3000};
    std::string serialized_data = oss.str();

    // Now, 'serialized_data' contains the serialized representation of the PrimalGrid
    // You can send 'serialized_data' using MPI or store it in a file, etc.

    // For deserialization, you can reconstruct the PrimalGrid from the bytes as follows:
    PrimalGrid grid_deserialized;
    std::istringstream iss(serialized_data);
    boost::archive::binary_iarchive ia(iss);
    ia >> grid_deserialized;

    // Now 'grid_deserialized' contains the deserialized PrimalGrid with the original data

    cout << "second_print\n";
    grid_deserialized.print();

    return 0;
}
