#pragma once
#include "../geometry/FV_Grid.hpp"
#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

/*--------------------------------------------------------------------
Defining serialization for some classes in order to pack their content
into bytes and send to other processes with MPI. Most classes will
have the serialization function defined in the class header
--------------------------------------------------------------------*/

namespace serialization
{
    template <typename T>
    inline void serialize(string &bytes, const T &obj)
    {
        ostringstream oss;
        boost::archive::binary_oarchive oa(oss);
        oa << obj;
        bytes = oss.str();
    }
    template <typename T>
    inline void deserialize(const string &input_bytes, T &obj)
    {
        istringstream iss(input_bytes);
        boost::archive::binary_iarchive ia(iss);
        ia >> obj;
    }
}

namespace boost::serialization
{
    template <class Archive>
    void serialize(Archive &ar, Vec3 &v, const unsigned int version)
    {
        ar &v[0];
        ar &v[1];
        ar &v[2];
    }

    // template <class Archive>
    // void serialize(Archive &ar, geometry::ElementPatch &obj, const unsigned int version) {
    // }

    // template <class Archive>
    // void serialize(Archive &ar, geometry::PatchBoundary &obj, const unsigned int version) {}

    // template <class Archive>
    // void serialize(Archive &ar, geometry::PatchInterface &obj, const unsigned int version) {}

}