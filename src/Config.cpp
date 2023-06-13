
#include "../include/Config.hpp"

Config::Config(string config_filename)      
{
    //Todo: fix config file
    

    //tmp
    mesh_filename = "/home/anders/dev/Compress3D/meshing/brick.su2";
    map_patch_BC = {{"inlet", BoundaryType::NoSlipWall},
                        {"outlet", BoundaryType::SlipWall},
                         {"sides", BoundaryType::NoSlipWall}};
    
    if (time_scheme == TimeScheme::ExplicitEuler || time_scheme == TimeScheme::TVD_RK3)
        time_integration_type = TimeIntegrationType::Explicit;
    else{
        FAIL_MSG("Implicit schemes are not yet implemented\n");
    }

}

void Config::set_grid_metrics(Index N_NODES, 
                    Index N_INTERIOR_CELLS, 
                    Index N_TOTAL_CELLS, 
                    Index N_INTERIOR_FACES,
                    Index N_TOTAL_FACES){
    assert(!grid_metrics_set);
    grid_metrics_set = true;

    this->N_NODES = N_NODES;
    this->N_TETS = N_INTERIOR_CELLS;
    this->N_INTERIOR_CELLS = N_INTERIOR_CELLS;
    this->N_TOTAL_CELLS = N_TOTAL_CELLS;
    this->N_INTERIOR_FACES = N_INTERIOR_FACES;
    this->N_TOTAL_FACES = N_TOTAL_FACES;
}


ConfigParser::ConfigParser(string config_filename){
    try {
        cfg_node = YAML::LoadFile(config_filename);
    } catch (YAML::Exception& e){
        FAIL_MSG("Error loading config file " << e.what());
    }
}


void ConfigParser::read_options(){

    for (const auto& pair : cfg_node){
        string key = pair.first.as<string>();
        std::type_index value_type = data_type_map.at(key);

        

    }
}