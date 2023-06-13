
#include "test.hpp"
#include <bit>

#include <yaml-cpp/yaml.h>

#include <exception>

struct Config{

    string solverName;
    int maxIterations;
    int tolerance;
    string algo;



    void load(string filePath){

        std::ifstream file(filePath);
        std::stringstream buffer;
        buffer << file.rdbuf();
        file.close();

        YAML::Node config = YAML::Load(buffer.str());
        solverName = config["solver"]["name"].as<std::string>();
        maxIterations = config["solver"]["max_iterations"].as<int>();
        tolerance = config["solver"]["tolerance"].as<double>();
        algo = config["algo"]["name"].as<string>();


    }



    void try_load(string path){
        try{
            load(path);
        }catch (YAML::Exception& e){
            std::cerr << "file read failed, " << e.msg << std::endl;
            exit(1);
        }
    }




};

void what_is_x(int x){

    try{
    
    if (x==0)
        throw std::runtime_error("error: x == 0");
    assert(x != 1);

    } catch(const std::exception& e){
        std::cerr << e.what() <<endl;
    }            
}


int main()
{

    //what_is_x(1);

    std::vector v = {1,2,3};
    cout << v.at(4);
}