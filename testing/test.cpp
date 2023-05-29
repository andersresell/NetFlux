
#include "test.hpp"

using namespace std;
    using EMat = Eigen::Matrix<double, 4,1>;

// inline void convert_to_static_type2(Index i, Static& s, double* data)
// {
//     //static_assert(sizeof(StaticContainerType) == 4*sizeof(double));


//     //s = *(Static*)&data[i*4];

//     s = *reinterpret_cast<Static*>(&data[i * 4]);
// }
template<typename T>
void print(T* d, int n){
    for (int i{0};i<n;i++) cout << d[i] <<", ";
    cout << endl<<endl;
}

void foo(Eigen::Map<EMat>& mat){
    mat *= 3;
    cout << "mat:" <<mat<<endl<<endl;
}

int main()
{
    
    Dynamic d;




    Eigen::Map<EMat> mat = d.get_variable<EMat>(1);
    foo(mat);
    cout << mat << endl;
    mat[0] = 50;

    d.print();
    
    mat = d.get_variable<EMat>(0);

    cout << mat<<endl;
    // EMat mat = d.get_variable<EMat>(1);
    // cout << mat << endl;

    // mat *= 0;

    // d.print();

    // int N{16};
    
    // double* data= new double[N];
    // for (int i{0};i<N;i++) data[i] = i;

    // print(data, N);

    // using EMat = Eigen::Matrix<double, 4,1>;

    // //EMat mat;
    // //EMat& mat = Eigen::Map<EMat&>(data);

    // Eigen::Matrix<double,4,1> mat1;
    // mat1.setConstant(1);

    // Eigen::Matrix<double, 4,1>::Map(data) = mat1;

    // data[0] = 0;
    // cout << mat1<<"\n\n";


    // print(data,N);
    using namespace Eigen;

    int data[] = {1,2,3,4,5,6,7,8,9};
    Map<RowVectorXi> v(data,4);
    cout << "The mapped vector v is: " << v << "\n";
    new (&v) Map<RowVectorXi>(data+4,5);
    cout << "Now v is: " << v << "\n";
    print<int>(data,9);

    // using EMat = Eigen::Matrix<int, 4,1>;
    // Map<EMat> mat(data);
    // cout << "mat "<<"\n"<<mat<<endl;
    // mat[0] = 200;
    // print<int>(data,9);


}