#pragma once
#include "Config.hpp"
#include <typeindex>
#include <typeinfo>

class ConfigParser{

    YAML::Node cfg_node;


public:
    ConfigParser(string config_filename);

    /*Sets member variables in the Config object*/
    
    void parse_config(Config& config);


private:

    // const std::map<string, std::type_index> data_type_map{
    //     {"mesh_filename", std::type_index(typeid(string))},
    //     {"output_basename", std::type_index(typeid(string))},
    //     {"solver", std::type_index(typeid(int))}
    // };

    /*Specifies options that are not specified in the 
    config file, but can be inferred after config file options are set*/
    void infer_hidden_options(Config& config);

    template<typename T>
    T read_required_option(string option_name){
        if (cfg_node[option_name]){
            return cfg_node[option_name].as<T>(); 
        } else {
            throw std::runtime_error(option_name + " not specified");
        }
    }

    template<typename T>
    T read_optional_option(string option_name, T default_value){
        if (cfg_node[option_name]){
            return cfg_node[option_name].as<T>(); 
        } else {
            return default_value;
        }
    }

    template<typename EnumType>
    EnumType read_required_enum_option(string option_name, const map<string, EnumType>& enum_map){
        if (cfg_node[option_name]){
            return enum_map.at(cfg_node[option_name].as<string>()); 
        } else {
            throw std::runtime_error(option_name + " not specified");
        }
    }

    template<typename EnumType>
    EnumType read_optional_enum_option(string option_name, const map<string, EnumType>& enum_map, EnumType default_value){
        if (cfg_node[option_name]){
            return enum_map.at(cfg_node[option_name].as<string>()); 
        } else {
            return default_value;
        }
    }

};