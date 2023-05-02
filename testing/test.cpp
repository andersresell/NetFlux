
#include "test.hpp"
#include <string>


int main(){
    string solver = "Euler";

    if (solver == "Euler"){
        EulerDriver ed();
        ed.solve();
    }


}