#include <iostream>
#include <vector>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <mpi.h>
#include <sstream>

// Define your data structure (example: MyClass)
class MyClass
{
public:
    int x;
    double y;

    // Add serialization function for MyClass (you can use Boost Serialization)
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &x;
        ar &y;
    }
};

using namespace std;

int main(int argc, char **argv)
{

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "hey rank " << rank << endl;

    MyClass my_obj;
    if (rank == 0)
    {
        vector<MyClass> objects;
        for (int i{0}; i < size; i++)
            objects.push_back({i + 1, 0.25 * (i + 1)});
        my_obj = move(objects[0]);

        for (int i{1}; i < size; i++)
        {
            ostringstream oss;
            boost::archive::binary_oarchive oa(oss);
            oa << objects[i];
            int num_bytes = oss.str().size();
            MPI_Send(&num_bytes, 1, MPI_INT, i, i, MPI_COMM_WORLD);
            MPI_Send(oss.str().data(), num_bytes, MPI_BYTE, i, i, MPI_COMM_WORLD);
        }
    }
    else
    {
        int num_bytes{-1};
        MPI_Recv(&num_bytes, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        string bytes_received;
        bytes_received.resize(num_bytes);
        MPI_Recv(bytes_received.data(), num_bytes, MPI_BYTE, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        istringstream iss(bytes_received);
        boost::archive::binary_iarchive ia(iss);
        ia >> my_obj;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "Finito, rank = " << rank << ", x = " << my_obj.x << ", y = " << my_obj.y << endl;
    MPI_Finalize();
}
