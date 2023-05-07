
#include "test.hpp"
#include <string>

#include <vector>
#include <iostream>
using namespace std;
int main(){
    A a;
    cout << a.a <<endl
         << a.b << endl
         << a.d <<endl
         << a.str <<endl;



    MyClass obj1{};
    
    std::cout << "Value 1: " << obj1.value << std::endl;
    cout << obj1.a.a <<endl
        << obj1.a.b << endl
        << obj1.a.d <<endl
        << obj1.a.str <<endl;
    return 0;
}