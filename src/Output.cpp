
#include "../include/Output.hpp"

Output::Output(const geom::Grid &grid, const Vector<unique_ptr<Solver>> &solvers, const Config &config)
    : grid{grid}, solvers{solvers}
{
    string output_dir = config.get_output_dir();

    if (filesys::exists(output_dir))
        if (!filesys::remove_all(output_dir))
            throw(std::runtime_error("Couldn't remove old output directory: " + output_dir));

    if (!filesys::create_directory(output_dir))
        throw std::runtime_error("Couldn't create output directory: " + output_dir);
}

void Output::write_vtk_ascii(const Config &config)
{
    /*Only write output every WRITE_STRIDE times*/
    if (config.get_timestep() % config.get_write_stride() != 0)
        return;

    const string &filename = config.get_unsteady_vtk_filename();

    write_vtk_ascii_grid(config, filename);

    for (const auto &solver : solvers)
    {

        const VecField &consvars = solver->get_solver_data().get_solution();

        switch (solver->get_solver_type())
        {
        case SolverType::Euler:
        {
            EulerOutput::write_vtk_ascii_cell_data(config, filename, consvars);
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
    std::ofstream ost{filename};
    FAIL_IF_MSG(!ost, "Couldn't open file " + filename);

    const Index N_NODES = config.get_N_NODES();
    const Index N_TETS = config.get_N_TETS();
    const ShortIndex VTK_TET_TYPE = 10; // VTK_TRI_TYPE = 5;
    const auto &nodes = grid.get_nodes();
    const auto &tet_connectivity = grid.get_tet_connect();

    ost << "# vtk DataFile Version 3.0\n"
        << "Compress 3D\n"
        << "ASCII\n\n"
        << "DATASET UNSTRUCTURED_GRID\n\n"
        << "POINTS " << N_NODES << " " + string(Scalar_name) + "\n";

    /*Writing the grid*/
    for (const auto &node : nodes)
        ost << node.x() << " " << node.y() << " " << node.z() << "\n";

    ost << "\nCELLS " << N_TETS << " " << (1 + N_TET_FACES) * N_TETS << "\n";
    for (const auto &tc : tet_connectivity)
        ost << N_TET_FACES << " " << tc.a() << " " << tc.b() << " " << tc.c() << " " << tc.d() << "\n";

    ost << "\nCELL_TYPES " << N_TETS << "\n";
    for (Index i{0}; i < N_TETS; i++)
        ost << VTK_TET_TYPE << "\n";
}

/*Writing various fields such as limiters or gradients*/
void Output::write_vtk_ascii_debug(const Config &config, const string &filename)
{
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

void EulerOutput::write_vtk_ascii_cell_data(const Config &config, const string &filename, const VecField &consvars)
{
    assert(consvars.get_N_EQS() == N_EQS_EULER);

    std::ofstream ost{filename, std::ios::app};

    if (!ost)
        throw std::runtime_error("Couldn't open file " + filename);

    const Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();

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