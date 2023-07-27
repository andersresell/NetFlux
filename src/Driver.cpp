#include "../include/Driver.hpp"
#include "../include/parallelization/DomainDecomposition.hpp"

Driver::Driver(Config &config) : config{config}
{
    cout << "\n///////////////////////////////////\n"
         << "////////  ->| NetFlux |->  ////////\n"
         << "///////////////////////////////////\n\n";

    create_partitioned_grids(config, primal_grid, FV_grid);

    switch (config.get_main_solver_type())
    {
    case MainSolverType::Euler:
        solvers.push_back(make_unique<EulerSolver>(config, *primal_grid, *FV_grid));
        break;
    default:
        throw std::runtime_error("Error: Illegal solver type specified");
    }

    output = make_unique<Output>(*primal_grid, solvers, config);
}

void Driver::solve()
{
    cout << "\n----- Starting solver -----\n";
    config.start_counter();
    config.set_timestep(0);
    config.set_time(0.0);
    output->write_vtk_ascii(config);
    while (1)
    {
        /*If main solvers consisting of multiple sub-solvers (for instance NS + scalar transport) are implemented in the future,
        additional constructs should be added to handle the sequential solution of those. This is ignored for now. It is
        assumed that there is only one solver*/
        assert(solvers.size() == 1);

        solvers.at(0)->calc_timestep(config);

        for (auto &solver : solvers)
        {
            solver->step(config);
        }

        config.set_timestep(config.get_timestep() + 1);

        config.set_time(config.get_time() + config.get_delta_time());

        cout << "Time step " << config.get_timestep() << " finished\n";

        output->write_vtk_ascii(config);

        if (config.get_timestep() >= config.get_n_timesteps())
        {
            break;
        }
    }

    cout << "Solver finished\n"
         << "Elapsed time: " << config.get_elapsed_time() << endl;
}
