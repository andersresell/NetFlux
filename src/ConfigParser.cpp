#include "../include/ConfigParser.hpp"

ConfigParser::ConfigParser(string config_filename){
    try {
        cfg_node = YAML::LoadFile(config_filename);
    } catch (YAML::Exception& e){
        FAIL_MSG("Error loading config file " << e.what());
    }
}

    
void ConfigParser::parse_config(Config& config){
    
    
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
    
    
}

void ConfigParser::infer_hidden_options(Config& config){
    
    if (config.time_scheme == TimeScheme::ExplicitEuler || config.time_scheme == TimeScheme::TVD_RK3){
        config.time_integration_type = TimeIntegrationType::Explicit;
    
    }else{
        FAIL_MSG("Implicit schemes are not yet implemented\n");
    }

    
}