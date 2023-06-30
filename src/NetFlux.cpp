

#include "../include/Driver.hpp"

int main(int argc, char *argv[])
{
    try
    {
        if (argc != 2)
        {
            cout << "Specify config file to run NetFlux.\nUsage: $ ./NetFlux <config_file>\n";
            return 1;
        }
        string config_filename = argv[1];
        Config config{config_filename};
        Driver driver{config};
        driver.solve();
    }
    catch (std::exception &e)
    {
        std::cerr << "Exception caught:\n"
                  << e.what() << endl;
        return 1;
    }

    return 0;
}