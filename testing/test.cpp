
#include "test.hpp"
#include <bit>

using namespace std;
    using EMat = Eigen::Matrix<double, 4,1>;
    using EVec = Eigen::Vector<double, 4>;
// 

struct FlowVar{
    double density,u,v,w,p;
    void print(){
        cout << density<<", "<<u<<", "<<v<<", "<<w<<", "<< p<<endl;
    }
};

    template<size_t N_VARS>
    struct FlowVars{
        double variables [N_VARS];
    };

int main()
{

    double arr [10];
    for (Index i{0};i<10;i++)arr[i]=i;

    const double* ptr = arr;

    for (Index i{0};i<13;i++)cout<<ptr[i]<<endl;

    EMat e{ptr};
    cout << e;
    arr[1] = 121221;
    cout << e << endl;

    Eigen::Map<EMat> emap{e.data()};

    cout << sizeof(emap)<<endl;


    cout << sizeof(Eigen::Map<EVec>)<<endl;

    
}