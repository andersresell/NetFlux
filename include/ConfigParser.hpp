#pragma once
#include "Config.hpp"

namespace YAML
{

    template <>
    struct convert<Vec3>
    {
        static bool decode(const Node &node, Vec3 &vec)
        {
            if (!node.IsSequence())
            {
                return false;
            }

            if (node.size() != Vec3::RowsAtCompileTime)
                throw std::runtime_error("The vector " + node.Tag() + " has the wrong length\n");

            size_t i{0};
            for (const auto &element : node)
            {
                vec[i++] = element.as<Scalar>();
            }

            return true;
        }
    };

}

class ConfigParser
{

    YAML::Node root_node;

    string config_filename;

    string sim_dir_path;

public:
    ConfigParser(const string &sim_dir);

    /*Sets member variables in the Config object*/

    void parse_config(Config &config);

private:
    // const std::map<string, std::type_index> data_type_map{
    //     {"mesh_filename", std::type_index(typeid(string))},
    //     {"output_basename", std::type_index(typeid(string))},
    //     {"solver", std::type_index(typeid(int))}
    // };

    /*Specifies options that are not specified in the
    config file, but can be inferred after config file options are set*/

    void parse_yaml_file_options(Config &config);

    void infer_hidden_options(Config &config);

    string option_not_specified_msg(string option_name) const
    {
        return "\"" + option_name + "\" not specified in the configuration file \"" + config_filename + "\"";
    }

    template <typename EnumType>
    EnumType lookup_enum_option_map(const map<string, EnumType> &map,
                                    const string &key,
                                    string option_name,
                                    const string &option_parent_name = "") const
    {
        if (map.count(key) == 1)
            return map.at(key);
        else
        {
            string keys;
            for (const auto &pair : map)
                keys += "'" + pair.first + "'\n";
            if (option_parent_name.size() > 0)
                option_name = option_parent_name + ": " + option_name;
            throw std::runtime_error("Illegal value '" + key + "' specified for setting '" +
                                     option_name + "'. Legal values are:\n" + keys);
        }
    }

    template <typename T>
    T read_required_option(string option_name)
    {
        if (root_node[option_name])
        {
            return root_node[option_name].as<T>();
        }
        else
        {
            throw std::runtime_error(option_not_specified_msg(option_name));
        }
    }

    template <typename T>
    T read_optional_option(string option_name, T default_value)
    {
        if (root_node[option_name])
        {
            return root_node[option_name].as<T>();
        }
        else
        {
            return default_value;
        }
    }

    template <typename EnumType>
    EnumType read_required_enum_option(string option_name, const map<string, EnumType> &enum_map)
    {
        if (root_node[option_name])
        {
            return lookup_enum_option_map(enum_map, root_node[option_name].as<string>(), option_name);
        }
        else
        {
            throw std::runtime_error(option_not_specified_msg(option_name));
        }
    }

    template <typename EnumType>
    EnumType read_optional_enum_option(string option_name, const map<string, EnumType> &enum_map, EnumType default_value)
    {
        if (root_node[option_name])
        {
            return lookup_enum_option_map(enum_map, root_node[option_name].as<string>(), option_name);
        }
        else
        {
            return default_value;
        }
    }

    void read_PatchExtes(Config &config);

    void set_mesh_name(Config &config);
};