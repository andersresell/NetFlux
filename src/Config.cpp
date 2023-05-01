
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
        FAIL_MSG("Implicit schemes not yet implemented\n");
    }

}

void Config::set_grid_data(Index N_NODES, 
                    Index N_INTERIOR_CELLS, 
                    Index N_TOTAL_CELLS, 
                    Index N_FACES){
    assert(!grid_data_set);
    grid_data_set = true;

    this->N_NODES = N_NODES;
    this->N_TETS = N_INTERIOR_CELLS;
    this->N_INTERIOR_CELLS = N_INTERIOR_CELLS;
    this->N_TOTAL_CELLS = N_TOTAL_CELLS;
    this->N_FACES = N_FACES;
}