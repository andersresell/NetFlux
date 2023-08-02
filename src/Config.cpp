
#include "../include/Config.hpp"
#include "../include/ConfigParser.hpp"

Config::Config(string sim_dir_path)
{
    ConfigParser parser{sim_dir_path};
    parser.parse_config(*this);
}

void Config::set_grid_metrics_local(Index N_NODES,
                                    Index N_INTERIOR_CELLS,
                                    Index N_GHOST_CELLS_PART,
                                    Index N_GHOST_CELLS_EXT,
                                    Index N_INTERIOR_FACES,
                                    Index N_PARTITION_FACES,
                                    Index N_EXTERIOR_FACES)
{
    assert(!grid_metrics_loc_set);
    this->N_NODES = N_NODES;
    this->N_CELLS_INT = N_INTERIOR_CELLS;
    this->N_CELLS_GHOST_PART = N_GHOST_CELLS_PART;
    this->N_CELLS_GHOST_EXT = N_GHOST_CELLS_EXT;
    this->N_FACES_INT = N_INTERIOR_FACES;
    this->N_FACES_PART = N_PARTITION_FACES;
    this->N_FACES_EXT = N_EXTERIOR_FACES;
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
    this->N_CELLS_INT_GLOB = N_INTERIOR_CELLS;
    this->N_CELLS_GHOST_GLOB = N_TOTAL_CELLS - N_INTERIOR_CELLS;
    this->N_FACES_INT_GLOB = N_INTERIOR_FACES;
    this->N_FACES_EXT_GLOB = N_TOTAL_FACES - N_INTERIOR_FACES;
    grid_metrics_glob_set = true;
}

bool Config::valid_mesh_name(const string &name) const
{
    for (const auto &extension : valid_mesh_extensions)
        if ((string(MESH_NAME_NO_EXTENSION) + "." + string(extension)) == name)
            return true;
    return false;
}

BoundaryType Config::get_boundary_type(const string &PatchExt_name) const
{
    assert(map_PatchExt_BC.count(PatchExt_name) == 1); // This should have been caught when reading mesh file
    return map_PatchExt_BC.at(PatchExt_name);
}