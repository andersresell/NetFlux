
#include "../include/Config.hpp"
#include "../include/ConfigParser.hpp"

Config::Config(string sim_dir_path)
{
    ConfigParser parser{sim_dir_path};
    parser.parse_config(*this);
}

void Config::set_grid_metrics(Index N_NODES,
                              Index N_INTERIOR_CELLS,
                              Index N_TOTAL_CELLS,
                              Index N_INTERIOR_FACES,
                              Index N_TOTAL_FACES)
{
    assert(!grid_metrics_set);
    grid_metrics_set = true;

    this->N_NODES = N_NODES;
    this->N_TETS = N_INTERIOR_CELLS;
    this->N_INTERIOR_CELLS = N_INTERIOR_CELLS;
    this->N_TOTAL_CELLS = N_TOTAL_CELLS;
    this->N_INTERIOR_FACES = N_INTERIOR_FACES;
    this->N_TOTAL_FACES = N_TOTAL_FACES;
}
