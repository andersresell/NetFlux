
#include "test.hpp"

using namespace std;
int main()
{
    

    using CMat = Container2D<double, 3, 5>;
    using EMat = Eigen::Matrix<double, 3, 5>;

    const Index M = 5, N = 3, K=1;

    Container2D<double,5,3> grad;
    grad=2;
    grad(4,0) = 1;
    using Vec3 = Eigen::Vector3d;
    Vec3 x;

    x[0] = 1;x[1] = 2; x[2] = 3;
    Container1D<double,5> U;

    container::Vec3_mult(grad,x,U);

    cout << "U:\n"<<U<<endl;


    
}