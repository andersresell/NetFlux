
#include "test.hpp"

template <size_t N,typename T>
class MyClass {
public:
    T data[N];

public:
    template <typename... Args>
    MyClass(Args... args) : data{static_cast<T>(args)...} {
        static_assert(sizeof...(args) == N);
    }

    // other member functions here
};

template <typename T>
class FixedSizeMyClass : public MyClass<2, T> {
public:
    template <typename... Args>
    FixedSizeMyClass(Args&&... args) : MyClass<2, T>(std::forward<Args>(args)...) {}
};



int main(){
    MyClass<2, double> obj0{1,2.1};

    for (int i{0}; i<2;i++) cout << obj0.data[i] << endl;


    Container<double,2> obj1{1,2.1};

    for (int i{0}; i<2;i++) cout << obj1[i] << endl;

    FixedSizeMyClass<double> fsm{1.0,3};
    for (int i{0}; i<2;i++) cout << fsm.data[i] << endl;

    EulerVar ev{2,3,4,3};
    cout << ev << endl;
}