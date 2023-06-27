// #include <iostream>
#include <cassert>

using namespace std;

int func1(int a)
{
    int b = a * a * a;
    assert(b > 0);
    return b;
}

int func2(int a)
{
    return a * a * a;
}
