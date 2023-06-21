

#include "../include/Driver.hpp"

int main(){


    // string mesh_filename = "/home/anders/dev/Compress3D/meshing/brick.su2";

    // Config config{mesh_filename};
    
    // geom::Grid g{c};  

    // g.create_grid(c);

    // g.print_grid(c);

    // Output o{g};
    try{

        string config_filename = "../testing/test.yaml";
        Config config{config_filename};
        Driver driver{config};
        driver.solve();
    } catch (std::exception& e)
    {
        std::cerr << "Exception caught:\n" << e.what() <<endl;
        return 1;
    }
    return 0;
}