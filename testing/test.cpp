
#include "test.hpp"

using namespace std;


inline void convert_to_static_type2(Index i, Static& s, double* data)
{
    //static_assert(sizeof(StaticContainerType) == 4*sizeof(double));


    //s = *(Static*)&data[i*4];

    s = *reinterpret_cast<Static*>(&data[i * 4]);
}


int main()
{
    
    Dynamic d;




    double* arr = new double[16];
    for (int i{0};i<16;i++) arr[i] = i;

    Static& s2 = *(Static*)&d.ddata[0];

    s2.print();
    for (int i{0};i<4;i++) (s2)[i] = 0;
    s2.print();

    d.print();

    // Static s1; 

    // d.print();

    // convert_to_static_type2(0,s1, d.ddata);
    // s1.print();
    // for (int i{0};i<4;i++) s1[i] = 0;
    // s1.print();
    // d.print();


}