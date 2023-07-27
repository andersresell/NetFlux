#include "../../include/parallelization/DomainDecomposition.hpp"

/*--------------------------------------------------------------------
This function reads the mesh file and creates both a primal grid
and a finite volume grid. If multiple processors are used, Metis
is used to partition the mesh.
--------------------------------------------------------------------*/

void create_partitioned_grids(Config &config,
                              unique_ptr<geometry::PrimalGrid> &primal_grid,
                              unique_ptr<geometry::FV_Grid> &FV_grid)
{
    using namespace geometry;
    const ShortIndex num_procs = NF_MPI::get_size();
    const ShortIndex rank = NF_MPI::get_rank();

    /*--------------------------------------------------------------------
    Rank 0 creates a global grid.
    --------------------------------------------------------------------*/
    unique_ptr<PrimalGrid> primal_grid_glob;
    if (rank == 0)
    {
        primal_grid_glob = make_unique<PrimalGrid>(config);

        /*--------------------------------------------------------------------
        If there is only one processor, assign the global grid as the only
        primal grid, construct finite volume grid, and return.
        --------------------------------------------------------------------*/
        if (num_procs == 1)
        {
            primal_grid = move(primal_grid_glob);
            FV_grid = make_unique<geometry::FV_Grid>(config, *primal_grid);
            return;
        }
    }

    /*--------------------------------------------------------------------
    If multiple processors are used, the domain must be partitioned into several
    pieces. This is handled by METIS.
    The return value 'volume_element_partition', is a vector of indices
    from 0 to n_procs-1 marking which processor each volume element
    belongs to.
    --------------------------------------------------------------------*/
    Vector<Index> volume_element_partition = NF_METIS::calc_element_partition(*primal_grid_glob, num_procs);

    const Elements &volume_elements_glob = primal_grid_glob->get_volume_elements();
    const Vector<Vec3> &nodes_glob = primal_grid_glob->get_nodes();
    assert(volume_element_partition.size() == volume_elements_glob.size());

    /*Count num nodes for each rank (only used for preallocation)*/
    Index num_volume_elements_loc{0};
    for (Index i : volume_element_partition)
        if (i == rank)
            num_volume_elements_loc++;

    Elements volume_elements_loc;
    volume_elements_loc.reserve(num_volume_elements_loc, MAX_NODES_VOLUME_ELEMENT);
    Vector<Vec3> nodes_loc;
    nodes_loc.reserve(1.5 * nodes_glob.size() / num_procs); /*Assuming equal partition times a factor as estimate for number of nodes for each processor*/

    map<Index, Index> glob_to_loc_node_id;

    /*--------------------------------------------------------------------
    Looping over volume elements to copy from global elements and nodes
    to local elements and nodes.
    --------------------------------------------------------------------*/
    for (Index i{0}; i < volume_elements_glob.size(); i++)
    {
        if (volume_element_partition[i] == rank)
        {
            const Index *element = volume_elements_glob.get_element_nodes(i);
            ShortIndex num_nodes = volume_elements_glob.get_n_element_nodes(i);
            ElementType e_type = volume_elements_glob.get_element_type(i);

            /*--------------------------------------------------------------------
            Loops over all global node indices of an element. If a new node index is
            discovered, the point corresponding to the node index is added and the
            a mapping from global indices are added.
            --------------------------------------------------------------------*/
            for (ShortIndex k{0}; k < num_nodes; k++)
            {
                Index node_id_glob = element[k];
                if (glob_to_loc_node_id.count(node_id_glob) == 0)
                {
                    nodes_loc.emplace_back(nodes_glob[node_id_glob]);
                    Index node_id_loc = nodes_loc.size();
                    glob_to_loc_node_id.emplace(node_id_glob, node_id_loc);
                }
            }
            /*--------------------------------------------------------------------
            Copies an element from the global to the local domain, and renumbers
            the node indices from global to local values.
            --------------------------------------------------------------------*/
            volume_elements_loc.add_element_local(e_type, element, glob_to_loc_node_id);
        }
    }

    /*--------------------------------------------------------------------
    Looping over patches and copying the patch elements and the patch name
    from global to local.
    --------------------------------------------------------------------*/
    Vector<ElementPatch> element_patches_loc;

    const Vector<ElementPatch> &element_patches_glob = primal_grid->get_element_patches();

    for (const auto &element_patch_glob : element_patches_glob)
    {
        const Elements &boundary_elements_glob = element_patch_glob.boundary_elements;
        const string &patch_name_glob = element_patch_glob.patch_name;
        element_patches_loc.emplace_back();
        element_patches_loc.back().patch_name = patch_name_glob;
        Elements &boundary_elements_loc = element_patches_loc.back().boundary_elements;
        boundary_elements_loc.reserve(boundary_elements_glob.size(), MAX_NODES_FACE_ELEMENT);

        /*--------------------------------------------------------------------
        Looping over all elements of a patch
        --------------------------------------------------------------------*/
        for (Index i{0}; i < boundary_elements_glob.size(); i++)
        {
            ElementType e_type = boundary_elements_glob.get_element_type(i);
            const Index *element = boundary_elements_glob.get_element_nodes(i);
            ShortIndex num_nodes = boundary_elements_glob.get_n_element_nodes(i);

            /*--------------------------------------------------------------------
            Using the global to local node index map (from the volume element loop)
            to check wether a boundary face belongs to the local domain. Keep in mind
            that this map only contains the nodes belonging to the local domain.
            If all global nodes are common, the global face has to be part of the
            local domain. In that case, the boundary element is added. Interfaces
            between local domains are not handled in this step.
            --------------------------------------------------------------------*/
            bool shared_face = true;
            for (ShortIndex k{0}; k < num_nodes; k++)
            {
                Index node_id_glob = element[k];
                if (glob_to_loc_node_id.count(node_id_glob) != 1)
                {
                    shared_face = false;
                    break;
                }
            }
            if (shared_face)
                boundary_elements_loc.add_element_local(e_type, element, glob_to_loc_node_id);
        }
        boundary_elements_loc.shrink_to_fit();
    }
    element_patches_loc.shrink_to_fit();
    volume_elements_loc.shrink_to_fit();
    nodes_loc.shrink_to_fit();

    primal_grid = make_unique<PrimalGrid>(nodes_loc, volume_elements_loc, element_patches_loc);
}

namespace NF_METIS
{
    Vector<Index> calc_element_partition(PrimalGrid &primal_grid,
                                         Index n_partitions)
    {
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
#ifndef NDEBUG
        options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO; // Enable debug info
#endif

        Elements &elements = primal_grid.get_volume_elements();

        Vector<Index> &e_ptr = elements.get_e_ptr();
        Vector<Index> &e_ind = elements.get_e_ind();

        idx_t n_elements = static_cast<idx_t>(elements.size());
        idx_t n_nodes = primal_grid.get_nodes().size();

        idx_t n_common = 3; // Minimum 3 common node to get a face/edge for a 3D mesh (A triangle)
        idx_t objval = -1;
        Vector<Index> element_partition(n_elements);
        Vector<Index> node_partition(n_nodes); // Don't really need this one, but I get segfault if I input a nullptr instead of this
        idx_t n_part_idx_t = static_cast<idx_t>(n_partitions);

        int status = METIS_PartMeshDual(&n_elements,
                                        &n_nodes,
                                        (idx_t *)e_ptr.data(),
                                        (idx_t *)e_ind.data(),
                                        nullptr,
                                        nullptr,
                                        &n_common,
                                        &n_part_idx_t,
                                        nullptr,
                                        nullptr,
                                        &objval,
                                        (idx_t *)element_partition.data(),
                                        (idx_t *)node_partition.data());

        assert(metis_statuses.count((rstatus_et)status) == 1);
        if (status != 1)
            throw std::runtime_error("Error partitioning mesh using METIS. METIS_PartMeshDual return message: " +
                                     metis_statuses.at(static_cast<rstatus_et>(status)) + "\n");
        return element_partition;
    }
}