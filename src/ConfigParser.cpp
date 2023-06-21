#include "../include/ConfigParser.hpp"


ConfigParser::ConfigParser(const string& config_filename) : config_filename{config_filename} {
    try {
        root_node = YAML::LoadFile(config_filename);

    } catch (const std::exception& e){
        throw std::runtime_error("Error loading config file" + string(e.what()));
    }
}

    
void ConfigParser::parse_config(Config& config){

    try {

        parse_yaml_file_options(config);

        infer_hidden_options(config);

    } 
    
    catch(std::exception& e){
        throw std::runtime_error("Error parsing config file:\n" + string(e.what()));

    } 

    cout << "Configuration file parsed without errors\n";
}


void ConfigParser::parse_yaml_file_options(Config& config){
    
    config.mesh_filename = read_required_option<string>("mesh_filename");

    config.output_basename = read_optional_option<string>("output_basename", "output");
    
    config.main_solver_type = read_required_enum_option<MainSolverType>("solver", main_solver_from_string);

    config.time_scheme = read_required_enum_option<TimeScheme>("time_scheme", time_scheme_from_string);

    config.inv_flux_scheme = read_required_enum_option<InviscidFluxScheme>("inviscid_flux_scheme", inviscid_flux_scheme_from_string);
    
    config.spatial_order = read_required_enum_option<SpatialOrder>("spatial_order", spatial_order_from_string);

    if (config.spatial_order == SpatialOrder::Second){
        config.limiter = read_required_enum_option<Limiter>("limiter", limiter_from_string);
    }

    config.grad_scheme = read_optional_enum_option<GradientScheme>("grad_scheme", gradient_scheme_from_string, GradientScheme::GreenGauss);
    
    config.initial_cond_option = read_required_enum_option<InitialConditionOption>("initial_cond", initial_condition_option_from_string);

    Scalar density_fs = read_optional_option<Scalar>("density_fs", standard_air::density);
    config.set_primvars_inf(primvars_index::Density, density_fs);
    
    Vec3 velocity_fs = read_optional_option<Vec3>("velocity_fs", Vec3{0.0, 0.0, 0.0});
    for (ShortIndex i{0}; i<N_DIM; i++) 
        config.set_primvars_inf(primvars_index::Velocity(i), velocity_fs[i]);
    
    Scalar pressure_fs = read_optional_option<Scalar>("pressure_fs", standard_air::pressure);
    config.set_primvars_inf(primvars_index::Pressure, pressure_fs);

    config.n_timesteps = read_required_option<size_t>("n_timesteps");

    config.CFL = read_required_option<Scalar>("CFL");

    read_patches(config);




}
void ConfigParser::infer_hidden_options(Config& config){
    
    if (config.time_scheme == TimeScheme::ExplicitEuler || config.time_scheme == TimeScheme::TVD_RK3){
        config.time_integration_type = TimeIntegrationType::Explicit;
    
    }else{
        assert_msg(false, "Implicit schemes are not yet implemented\n");
    }

    
}


void ConfigParser::read_patches(Config& config){
    if (root_node["patches"]){
        YAML::Node patches_node = root_node["patches"];

        if (!patches_node.IsSequence())
            throw std::runtime_error("\"patches\" setting in configuration file is formatted crrectly (not a valid sequence)");
        

        for (const auto& patch_it : patches_node){
            
            if (!patch_it.IsMap()) throw std::runtime_error("Entry in \"patches\" setting in configuration file is not a key-value pair");
            
            const auto& patch_kv = *patch_it.begin();

            string patch_name = patch_kv.first.as<string>();

            if (config.map_patch_BC.count(patch_name) > 0)
                throw std::runtime_error("Duplicate patch name \"" + patch_name + "\" specified in config file");

            string boundary_type_key = patch_kv.second.as<string>();

            BoundaryType bc_type = lookup_enum_option_map(boundary_type_from_string, boundary_type_key, patch_name, "patches");

            config.map_patch_BC.emplace(patch_name, bc_type);
        }
    
    }
    else{
        throw std::runtime_error("\"patches\" option needs to be included in the config file");
    }

    
    
}
