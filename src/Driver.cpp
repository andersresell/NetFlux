#include "../include/Driver.hpp"

Driver::Driver(Config& config){
    
    grid = std::make_unique<geom::Grid>(config);
    
    grid->create_grid(config);
    
    switch (config.get_main_solver_type()){
        case MainSolverType::Euler:
            solvers.push_back(make_unique<EulerSolver>(config, *grid));
            break;
        dafault:
            FAIL_MSG("Error: Illegal solver type specified");
    }

    output = make_unique<Output>(*grid, solvers);
        
}

void Driver::solve(Config& config){
    config.set_timestep(0);
    config.set_time(0.0);
    

    output->write_vtk_ascii(config);

    while (1){
        /*If main solvers consisting of multiple sub-solvers (for instance NS + scalar transport) are implemented in the future,
        additional constructs should be added to handle the sequential solution of those. This is ignored for now. It is
        assumed that there is only one solver*/
        assert(solvers.size() == 1);
        

        double delta_time = solvers[0]->calc_timestep(config);
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