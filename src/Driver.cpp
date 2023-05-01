#include "../include/Driver.hpp"

Driver::Driver(Config& config){
    
    grid = std::make_unique<geom::Grid>(config);
    
    grid->create_grid(config);
    
    solver = std::make_unique<EulerSolver>(config, *grid);
    
    output = std::make_unique<EulerOutput>(*grid, solver->get_solution());

    
}

void Driver::solve(Config& config){
    config.set_timestep(0);
    config.set_time(0.0);
    

    output->write_vtk_ascii(config);

    while (1){
        

        double delta_time = solver->calc_timestep(config);
        config.set_delta_time(delta_time);


        output->write_vtk_ascii(config);
    
        
        size_t& timestep = config.get_timestep();        
        timestep++;    

        double& time = config.get_time();
        time += delta_time;


        
        if (timestep > config.get_n_timesteps()){
            break;
        }

    }
    



    
}