
#include "../include/Config.hpp"
#include "../include/ConfigParser.hpp"

Config::Config(string sim_dir_path)
{
    ConfigParser parser{sim_dir_path};
    parser.parse_config(*this);
}

void Config::set_grid_metrics_local(Index N_NODES,
                                    Index N_INTERIOR_CELLS,
                                    Index N_TOTAL_CELLS,
                                    Index N_INTERIOR_FACES,
                                    Index N_TOTAL_FACES)
{
    assert(!grid_metrics_loc_set);
    this->N_NODES_LOC = N_NODES;
    this->N_INTERIOR_CELLS_LOC = N_INTERIOR_CELLS;
    this->N_TOTAL_CELLS_LOC = N_TOTAL_CELLS;
    this->N_INTERIOR_FACES_LOC = N_INTERIOR_FACES;
    this->N_TOTAL_FACES_LOC = N_TOTAL_FACES;
    grid_metrics_loc_set = true;
}
void Config::set_grid_metrics_global(Index N_NODES,
                                     Index N_INTERIOR_CELLS,
                                     Index N_TOTAL_CELLS,
                                     Index N_INTERIOR_FACES,
                                     Index N_TOTAL_FACES)
{
    assert(!grid_metrics_glob_set);
    this->N_NODES_GLOB = N_NODES;
    this->N_INTERIOR_CELLS_GLOB = N_INTERIOR_CELLS;
    this->N_TOTAL_CELLS_GLOB = N_TOTAL_CELLS;
    this->N_INTERIOR_FACES_GLOB = N_INTERIOR_FACES;
    this->N_TOTAL_FACES_GLOB = N_TOTAL_FACES;
    grid_metrics_glob_set = true;
}

bool Config::valid_mesh_name(const string &name) const
{
    for (const auto &extension : valid_mesh_extensions)
        if ((string(MESH_NAME_NO_EXTENSION) + "." + string(extension)) == name)
            return true;
    return false;
}

BoundaryType Config::get_boundary_type(const string &patch_name) const
{
    assert(map_patch_BC.count(patch_name) == 1); // This should have been caught when reading mesh file
    return map_patch_BC.at(patch_name);
}