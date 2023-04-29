
#include "test.hpp"

int main(){
    using namespace std;

    
    //Vector<int> v = {1,2,3};
    //v[3] = 1;

    Container<double,4> c;
    cout << "c "<< c << endl;

    Container<int,5> c1{1,2,3,4,5};
    Container<int,5> c2 = c1;
    c2[0] = 0;
    cout << c1<< endl;
    c1 *= 3;
    cout << c1 << endl;
    c1 -= c1;
    cout << c1 << endl;
    c1 -= c2;
    cout << c1 << endl;
    c2 = c1*4;
    
    Container<Index,4> rhs;
    TetConnect tc = {1,2,3,4};
    cout << tc << endl;
    tc.a() = 1;
    cout << tc << endl;
}