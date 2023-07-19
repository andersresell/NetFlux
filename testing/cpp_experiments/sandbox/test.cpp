#include <iostream>
#include <cassert>
#include <limits>
using namespace std;
#include "../../include/Includes.hpp"
#include "test.hpp"

#include <iostream>
#include <vector>
#include <algorithm>

struct MyType
{
    int value;
    // Add comparison operator for sorting
    bool operator<(const MyType &other) const
    {
        return value < other.value;
    }
};

struct S
{

    ShortIndex c;
    bool a, b, d;
};

int main()
{
    cout << "sz Index " << sizeof(Index) << endl;
    cout << "sz CellPair " << sizeof(Faces::CellPair) << endl;
    cout << "sz ShortIndex " << sizeof(ShortIndex) << endl;
    cout << "sz size_t " << sizeof(size_t) << endl;
    cout << "sz S " << sizeof(S) << endl;
    cout << "sz Vec3 " << sizeof(Vec3) << endl;
}
