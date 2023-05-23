#pragma once
#include "includes.hpp"
#include "Utilities.hpp"

class Config {

    Index N_NODES,
          N_TETS;

    Index N_INTERIOR_CELLS,
          N_TOTAL_CELLS,
          N_INTERIOR_FACES,
          N_TOTAL_FACES;
    bool grid_metrics_set{false};

    string mesh_filename; 
    string output_basename;

    map<string, BoundaryType> map_patch_BC; //map from each patch to the bc type applied

    MainSolverType main_solver_type;

    TimeIntegrationType time_integration_type;
    
    TimeScheme time_scheme;

    SpatialOrder spatial_order;

    GradientScheme grad_scheme;

    InviscidFluxScheme inv_flux_scheme;

    Limiter limiter;

    size_t n_timesteps;

    size_t timestep;

    double delta_time;

    double time;

    double CFL;


public:
    Config(string config_filename); 

    Index get_N_NODES() const {return N_NODES;}
    Index get_N_TETS() const {return N_TETS;}
    Index get_N_INTERIOR_CELLS() const {return N_INTERIOR_CELLS;}
    Index get_N_TOTAL_CELLS() const { return N_TOTAL_CELLS;}
    Index get_N_GHOST_CELLS() const { return N_TOTAL_CELLS - N_INTERIOR_CELLS;}
    Index get_N_INTERIOR_FACES() const {return N_INTERIOR_FACES;}
    Index get_N_TOTAL_FACES() const {return N_TOTAL_FACES;}
    Index get_N_BOUNDARY_FACES() const {return N_TOTAL_FACES - N_INTERIOR_FACES;}
    
    

    void set_grid_metrics(Index N_NODES, 
                    Index N_INTERIOR_CELLS, 
                    Index N_TOTAL_CELLS, 
                    Index N_INTERIOR_FACES,
                    Index N_TOTAL_FACES);
    
    string get_mesh_filename() const {return mesh_filename;}
    string get_unsteady_vtk_filename() const {return output_basename + std::to_string(timestep) + ".vtk";}

    BoundaryType get_boundary_type(string patch_name) const {return map_patch_BC.at(patch_name);}

    MainSolverType get_main_solver_type() const {return main_solver_type;}
    void set_main_solver_type(MainSolverType val) {main_solver_type = val;}

    TimeIntegrationType get_time_integration_type() const {return time_integration_type;}

    TimeScheme get_time_scheme() const {return time_scheme;}
    void set_time_scheme(TimeScheme val) {time_scheme = val;}

    GradientScheme get_grad_scheme() const {return grad_scheme;}
    void set_grad_scheme(GradientScheme val) {grad_scheme = val;}

    SpatialOrder get_spatial_order() const {return spatial_order;}
    void set_spatial_order(SpatialOrder val) {spatial_order = val;}

    InviscidFluxScheme get_inv_flux_scheme() const {return inv_flux_scheme;}
    void set_inv_flux_scheme(InviscidFluxScheme val) {inv_flux_scheme = val;}

    Limiter get_limiter() const {return limiter;}
    void set_limiter(Limiter val) {limiter = val;}

    size_t get_n_timesteps() const {return n_timesteps;}

    size_t& get_timestep() {return timestep;}
    void set_timestep(size_t val) {timestep = val;}

    double& get_time() {return time;}
    void set_time(double val) {time = val;}

    double get_delta_time() const {return delta_time;}
    void set_delta_time(double val) {delta_time = val;}

    double get_CFL() const {return CFL;}
    void set_CFL(double val) {CFL = val;}
};

