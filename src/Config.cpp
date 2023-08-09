
#include "../include/Config.hpp"
#include "../include/ConfigParser.hpp"

Config::Config(string sim_dir_path)
{
    ConfigParser parser{sim_dir_path};
    parser.parse_config(*this);
}

void Config::set_grid_metrics_local(Index N_NODES,
                                    Index N_CELLS_INT,
                                    Index N_FACES_INT,
                                    Index N_FACES_PART,
                                    Index N_FACES_EXT)
{
    assert(!grid_metrics_loc_set);
    this->N_NODES = N_NODES;
    this->N_CELLS_INT = N_CELLS_INT;
    this->N_FACES_INT = N_FACES_INT;
    this->N_FACES_PART = N_FACES_PART;
    this->N_FACES_EXT = N_FACES_EXT;
    grid_metrics_loc_set = true;
}
void Config::set_grid_metrics_global(Index N_NODES_GLOB,
                                     Index N_CELLS_INT_GLOB,
                                     Index N_FACES_INT_GLOB,
                                     Index N_FACES_EXT_GLOB,
                                     Index N_CONNECTIVITY_INDICES_GLOB)
{
    assert(!grid_metrics_glob_set);
    this->N_NODES_GLOB = N_NODES_GLOB;
    this->N_CELLS_INT_GLOB = N_CELLS_INT_GLOB;
    this->N_FACES_INT_GLOB = N_FACES_INT_GLOB;
    this->N_FACES_EXT_GLOB = N_FACES_EXT_GLOB;
    this->N_CONNECTIVITY_INDICES_GLOB = N_CONNECTIVITY_INDICES_GLOB;
    grid_metrics_glob_set = true;
}

bool Config::valid_mesh_name(const string &name) const
{
    for (const auto &extension : valid_mesh_extensions)
        if ((string(MESH_NAME_NO_EXTENSION) + "." + string(extension)) == name)
            return true;
    return false;
}

BoundaryType Config::get_boundary_type(const string &PatchBoundary_name) const
{
    assert(map_PatchBoundary_BC.count(PatchBoundary_name) == 1); // This should have been caught when reading mesh file
    return map_PatchBoundary_BC.at(PatchBoundary_name);
}

/*Making sure that all local partitions have the same global values.
Only doing this for delta_time now, but might add more parameters later*/
void Config::communicate_global_values()
{
    ShortIndex size = NF_MPI::get_size();
    ShortIndex rank = NF_MPI::get_rank();

    Scalar dt_loc_send = delta_time;
    vector<Scalar> dt_loc_recv;
    if (rank == 0)
        dt_loc_recv.resize(size);

    NF_MPI::Gather(&dt_loc_send, dt_loc_recv.data(), 1, 0);

    Scalar dt_glob{0.0};
    if (rank == 0)
    {
        for (Scalar dt_loc : dt_loc_recv)
        {
            assert(dt_loc > 0.0);
            dt_glob = max(dt_glob, dt_loc);
        }
    }
    NF_MPI::Bcast(&dt_glob, 1, 0);
    assert(dt_glob >= delta_time);
}
