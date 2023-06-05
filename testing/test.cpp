
#include "test.hpp"

using namespace std;
    using EMat = Eigen::Matrix<double, 4,1>;

class A{
    virtual void f() {cout << a<< "ho\n"; };
    int g() {return 2;}
protected:
    int a{1};
public:
    void call_it(){this->f();}

};

class B : A{
public:
    virtual void f() {cout <<g()<< a<< "hey\n"; };
};

int main()
{
    A* a = new B{};
    a->call_it();

}