#pragma once
#include "includes.hpp"
#include "Utilities.hpp"

class Config {

    Index N_NODES,
          N_TETS;

    Index N_INTERIOR_CELLS,
          N_TOTAL_CELLS,
          N_FACES;
    bool grid_data_set{false};

    string mesh_filename; 
    string output_basename;

    map<string, BoundaryType> map_patch_BC; //map from each patch to the bc type applied

    TimeScheme time_scheme;

    GradientScheme grad_scheme;

    ConvScheme conv_scheme;

    Limiter limiter;

    size_t n_timesteps;

    size_t timestep;

    double delta_time;


public:
    Config(string config_filename); 

    Index get_N_NODES() const {return N_NODES;}
    Index get_N_TETS() const {return N_TETS;}
    Index get_N_INTERIOR_CELLS() const {return N_INTERIOR_CELLS;}
    Index get_N_TOTAL_CELLS() const { return N_TOTAL_CELLS;}
    Index get_N_GHOST_CELLS() const { return get_N_TOTAL_CELLS() - get_N_INTERIOR_CELLS();}
    Index get_N_FACES() const {return N_FACES;}

    void set_grid_data(Index N_NODES, 
                       Index N_INTERIOR_CELLS, 
                       Index N_TOTAL_CELLS, 
                       Index N_FACES);
    
    string get_mesh_filename() const {return mesh_filename;}
    string get_unsteady_vtk_filename() const {return output_basename + std::to_string(timestep) + ".vtk";}

    BoundaryType get_boundary_type(string patch_name) const {return map_patch_BC.at(patch_name);}


    TimeScheme get_time_scheme() const {return time_scheme;}
    void set_time_scheme(TimeScheme val) {time_scheme = val;}

    GradientScheme get_grad_scheme() const {return grad_scheme;}
    void set_grad_scheme(GradientScheme val) {grad_scheme = val;}

    ConvScheme get_conv_scheme() const {return conv_scheme;}
    void set_conv_scheme(ConvScheme val) {conv_scheme = val;}

    Limiter get_limiter() const {return limiter;}
    void set_limiter(Limiter val) {limiter = val;}

    size_t get_n_timesteps() const {return n_timesteps;}

    size_t& get_timestep() {return timestep;}
    void set_timestep(size_t val) {timestep = val;}
    
    double get_delta_time() const {return delta_time;}
    void set_delta_time(double val) {delta_time = val;}

};

