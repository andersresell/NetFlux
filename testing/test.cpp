#include <iostream>
#include <cassert>
#include <limits>
using namespace std;

int main()
{
    double a = 15.0 / 0.0;
    cout << "a " << a << endl;
    double result = std::numeric_limits<double>::max() * 0.0001; // Evaluates to +inf
    cout << result << endl;
}