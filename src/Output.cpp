
#include "../include/Output.hpp"

Output::Output(const unique_ptr<geometry::PrimalGrid> &primal_grid_glob,
               const unique_ptr<geometry::PrimalGrid> &primal_grid,
               const vector<unique_ptr<Solver>> &solvers,
               const Config &config)
    : primal_grid_glob_{*primal_grid_glob}, primal_grid_{*primal_grid}, solvers{solvers}
{
    if (NF_MPI::get_rank() == 0)
    {
        string output_dir = config.get_output_dir();

        if (filesys::exists(output_dir))
            if (!filesys::remove_all(output_dir))
                throw(std::runtime_error("Couldn't remove old output directory: " + output_dir));

        if (!filesys::create_directory(output_dir))
            throw std::runtime_error("Couldn't create output directory: " + output_dir);

        assert(solvers.size() == 1); // fix if adding more solvers
        for (ShortIndex i{0}; i < solvers.size(); i++)
        {
            Index N_CELLS = primal_grid_glob_.get_vol_elements().size();
            Index N_EQS = solvers[i]->get_solver_data().get_N_EQS();
            consvars_glob.emplace_back(N_CELLS, N_EQS);
        }
    }
    else
    {
        assert(primal_grid_glob.get() == nullptr);
    }
}

void Output::write_vtk_ascii(const Config &config)
{
    /*Only write output every WRITE_STRIDE times*/
    if (config.get_timestep() % config.get_write_stride() != 0)
        return;

    const string &filename = config.get_unsteady_vtk_filename();

    write_vtk_ascii_grid(config, filename);

    for (ShortIndex i{0}; i < solvers.size(); i++)
    {
        gather_cell_data(config, i);
        // const VecField &consvars = solver->get_solver_data().get_solution();

        switch (solvers[i]->get_solver_type())
        {
        case SolverType::Euler:
        {
            EulerOutput::write_vtk_ascii_cell_data(config, filename, consvars_glob[i]);
            break;
        }
        default:
            assert(false); // illegal solver type
        }
    }

    if (config.write_vtk_debug())
        write_vtk_ascii_debug(config, filename);
}

void Output::write_vtk_ascii_grid(const Config &config, const string &filename)
{
    if (NF_MPI::get_rank() != 0)
        return;
    std::ofstream ost{filename};
    FAIL_IF_MSG(!ost, "Couldn't open file " + filename);

    const Index N_NODES = config.get_N_NODES_GLOB();
    const Index N_CELLS = config.get_N_CELLS_INT_GLOB();
    const auto &nodes = primal_grid_glob_.get_nodes();
    const auto &vol_elements = primal_grid_glob_.get_vol_elements();
    assert(N_CELLS == vol_elements.size());

    // const vector<Index> &n_ptr = vol_elements.get_n_ptr();
    // const vector<Index> &n_ind = vol_elements.get_n_ind();
    const Index SUM_NODES_OVER_ALL_ELEMENTS = vol_elements.sum_nodes_over_all_elements();

    ost << "# vtk DataFile Version 3.0\n"
        << "NetFlux\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n\n"
        << "POINTS " << N_NODES << " " + string(Scalar_name) + "\n";

    /*Writing the grid*/
    for (const auto &node : nodes)
        ost << node.x() << " " << node.y() << " " << node.z() << "\n";

    ost << "\nCELLS " << N_CELLS << " " << N_CELLS + SUM_NODES_OVER_ALL_ELEMENTS << "\n";
    for (Index i{0}; i < N_CELLS; i++)
    {
        Index n_element_nodes = vol_elements.get_n_element_nodes(i);
        const Index *element_nodes = vol_elements.get_element_nodes(i);
        ost << n_element_nodes << " ";
        for (ShortIndex k{0}; k < n_element_nodes; k++)
            ost << element_nodes[k] << " ";
        ost << "\n";
    }

    ost << "\nCELL_TYPES " << N_CELLS << "\n";
    for (Index i{0}; i < N_CELLS; i++)
        ost << static_cast<ShortIndex>(vol_elements.get_element_type(i)) << "\n";
}

void Output::gather_cell_data(const Config &config, ShortIndex i_solver)
{
    /*Use mpi gather to copy all local cell data to a global container*/
    const VecField &U_loc = solvers[i_solver]->get_solver_data().get_solution();
    const Scalar *sendbuf = U_loc.data();
    Index count = U_loc.size() * U_loc.get_N_EQS();
    Scalar *recvbuf = nullptr;
    if (NF_MPI::get_rank() == 0)
    {
        VecField &U_glob = consvars_glob[i_solver];
        recvbuf = U_glob.data();
        assert(U_glob.size() == config.get_N_CELLS_INT_GLOB());
    }

    NF_MPI::Gather(sendbuf, recvbuf, count, 0);
}

void EulerOutput::write_vtk_ascii_cell_data(const Config &config, const string &filename, const VecField &consvars)
{
    if (NF_MPI::get_rank() != 0)
        return;
    assert(consvars.get_N_EQS() == N_EQS_EULER);

    std::ofstream ost{filename, std::ios::app};

    if (!ost)
        throw std::runtime_error("Couldn't open file " + filename);

    const Index N_INTERIOR_CELLS = config.get_N_CELLS_INT_GLOB();

    ost << "\nCELL_DATA " << N_INTERIOR_CELLS << "\n";

    ost << "\nSCALARS density " + string(Scalar_name) + " 1\n"
        << "LOOKUP_TABLE default\n";
    for (Index i{0}; i < N_INTERIOR_CELLS; i++)
        ost << consvars(i, 0) << "\n";

    ost << "\nVECTORS velocity " + string(Scalar_name) + "\n";
    for (Index i{0}; i < N_INTERIOR_CELLS; i++)
    {
        ost << consvars(i, 1) / consvars(i, 0) << " " << consvars(i, 2) / consvars(i, 0) << " " << consvars(i, 3) / consvars(i, 0) << "\n";
    }

    ost << "\nSCALARS pressure " + string(Scalar_name) + " 1\n"
        << "LOOKUP_TABLE default\n";

    for (Index i{0}; i < N_INTERIOR_CELLS; i++)
    {
        ost << EulerEqs::pressure(consvars.get_variable<EulerVec>(i)) << "\n";
    }
}

/*Writing various fields such as limiters or gradients*/
void Output::write_vtk_ascii_debug(const Config &config, const string &filename)
{
    assert(false); /*Needs to be adapted for parallel implementation*/

    std::ofstream ost{filename, std::ios::app};
    if (!ost)
        throw std::runtime_error("Couldn't open file " + filename);

    assert(solvers.size() == 1); // edit if changed
    const auto &solver_data = solvers[0]->get_solver_data();

    const auto &limiter = solver_data.get_primvars_limiter();
    const auto &gradient = solver_data.get_primvars_gradient();
    assert(gradient.get_N_EQS() == limiter.get_N_EQS());
    const ShortIndex N_EQS = limiter.get_N_EQS();
    assert(gradient.size() == limiter.size());
    assert(gradient.cols() == N_DIM);
    const Index N_INTERIOR_CELLS = limiter.size();

    for (ShortIndex k{0}; k < N_EQS; k++)
    {
        ost << "\nSCALARS limiter_EQ" + std::to_string(k) + " " + string(Scalar_name) + " 1\n"
            << "LOOKUP_TABLE default\n";
        for (Index i{0}; i < N_INTERIOR_CELLS; i++)
            ost << limiter(i, k) << "\n";
    }

    for (ShortIndex k{0}; k < N_EQS; k++)
    {
        ost << "\nVECTORS gradient_EQ" + std::to_string(k) + " " + string(Scalar_name) + "\n";
        for (Index i{0}; i < N_INTERIOR_CELLS; i++)
            ost << gradient(i, k, 0) << " " << gradient(i, k, 1) << " " << gradient(i, k, 2) << "\n";
    }
}
