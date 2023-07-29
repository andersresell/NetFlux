#include <iostream>
#include <vector>
#include <sstream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <eigen3/Eigen/Dense>

// Define a serialization function for Eigen::Vector3d
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

class PrimalGrid
{
public:
    std::vector<Eigen::Vector3d> nodes;

    // Define the serialization function for PrimalGrid
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &nodes;
    }
};
using namespace std;
int main()
{
    cout << "halla\n";
    // Create a PrimalGrid and add some Eigen::Vector3d objects
    PrimalGrid grid;
    grid.nodes.push_back(Eigen::Vector3d(1.0, 2.0, 3.0));
    grid.nodes.push_back(Eigen::Vector3d(4.0, 5.0, 6.0));
    grid.nodes.push_back(Eigen::Vector3d(7.0, 8.0, 9.0));

    // Serialize PrimalGrid to bytes
    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << grid;
    std::string serialized_data = oss.str();

    // Now, 'serialized_data' contains the serialized representation of the PrimalGrid
    // You can send 'serialized_data' using MPI or store it in a file, etc.

    // For deserialization, you can reconstruct the PrimalGrid from the bytes as follows:
    PrimalGrid grid_deserialized;
    std::istringstream iss(serialized_data);
    boost::archive::binary_iarchive ia(iss);
    ia >> grid_deserialized;

    // Now 'grid_deserialized' contains the deserialized PrimalGrid with the original data

    return 0;
}
