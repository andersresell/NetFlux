#include "../../include/parallelization/DomainDecomposition.hpp"

/*--------------------------------------------------------------------
This function reads the mesh file and creates both a primal grid
and a finite volume grid. If multiple processors are used, Metis
is used to partition the mesh.
--------------------------------------------------------------------*/

namespace geometry
{

    void GridCreator::construct_global_face_entities(const Elements &volume_elements_glob,
                                                     const Vector<ElementPatch> &element_patches,
                                                     map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                                     Elements &face_elements_glob,
                                                     map<FaceElement, InternalGhostData> &internal_boundary_faces,
                                                     const Vector<Index> &part)
    {
        if (NF_MPI::get_rank() != 0)
            return;

        assert(faces_to_cells_glob.size() == 0);

        for (Index cell_index{0}; cell_index < volume_elements_glob.size(); cell_index++)
        {
            ElementType volume_element_type = volume_elements_glob.get_element_type(cell_index);
            assert(is_volume_element(volume_element_type));
            for (ShortIndex k{0}; k < get_num_faces_volume_element(volume_element_type); k++)
            {
                FaceElement face_element = get_face_element_k_of_volume_element(volume_element_type,
                                                                                volume_elements_glob.get_element_nodes(cell_index),
                                                                                k);

                assert(faces_to_cells_glob.count(face_element) <= 1);
                if (faces_to_cells_glob.count(face_element) == 0)
                {
                    // Discovered a new face
                    faces_to_cells_glob.emplace(face_element, pair{cell_index, NOT_ASSIGNED});
                }
                else
                {
                    assert(faces_to_cells_glob.at(face_element).second == NOT_ASSIGNED);
                    // Face allready discovered by previous cell
                    faces_to_cells_glob.at(face_element).second = cell_index;
                }
            }
        }
        face_elements_glob.reserve(faces_to_cells_glob.size(), MAX_NODES_FACE_ELEMENT);

        /*--------------------------------------------------------------------
        Step 3: Adding internal faces to the map. (Boundary faces are added
        later, this is to get the correct grouping of patches). Face elements
        (used for calculating geometry properties) are created simultaneously
        as faces, to get the correct ordering.
        --------------------------------------------------------------------*/

        for (const auto &face : faces_to_cells_glob)
        {
            Index cell_i = face.second.first;
            long int cell_j = face.second.second;
            if (cell_j != NOT_ASSIGNED)
            { // Only add internal faces
                assert(cell_i < cell_j);
                face_elements_glob.add_element(face.first.e_type, face.first.nodes.data());

                /*Marking partition boundaries*/
                if (part[cell_i] != part[cell_j])
                {
                    assert(internal_boundary_faces.count(face.first) == 0);

                    internal_boundary_faces.emplace(face.first, InternalGhostData{cell_i,
                                                                                  (Index)cell_j,
                                                                                  (ShortIndex)part[cell_i],
                                                                                  (ShortIndex)part[cell_j]});
                }
            }
        }

        /*--------------------------------------------------------------------
        Step 4: Create the faces at the boundaries and ghost cells.
        Also add the face elements at the boundaries.
        --------------------------------------------------------------------*/

        // Index j_ghost = volume_elements_glob.size();
        for (const auto &element_patch : element_patches)
        {

            const Elements &surface_elements = element_patch.boundary_elements;
            for (Index ij{0}; ij < surface_elements.size(); ij++)
            {
                ElementType e_type = surface_elements.get_element_type(ij);
                FaceElement face_element{e_type, surface_elements.get_element_nodes(ij)};
                assert(faces_to_cells_glob.at(face_element).second == NOT_ASSIGNED);
                // faces_to_cells.at(face_element).second = j_ghost;
                Index i_domain = faces_to_cells_glob.at(face_element).first;
                face_elements_glob.add_element(e_type, face_element.nodes.data());
                // j_ghost++;
            }
        }
        assert(face_elements_glob.size() == faces_to_cells_glob.size());
    }

    void GridCreator::create_local_volume_entities(const Elements &volume_elements_glob,
                                                   const Vector<Vec3> &nodes_glob,
                                                   const Vector<ElementPatch> &element_patches_glob,
                                                   Index i_part,
                                                   const Vector<Index> &part,
                                                   Elements &volume_elements_loc,
                                                   Vector<Vec3> &nodes_loc,
                                                   Vector<ElementPatch> &element_patches_loc,
                                                   map<Index, Index> &nID_glob_to_loc,
                                                   map<Index, Index> &nID_loc_to_glob,
                                                   map<Index, Index> &eID_glob_to_loc)
    {

        if (NF_MPI::get_rank() != 0)
            return;

        assert(part.size() == volume_elements_glob.size());
        assert(volume_elements_loc.size() == 0);
        assert(nodes_loc.size() == 0);
        assert(element_patches_loc.size() == 0);
        assert(nID_glob_to_loc.size() == 0);

        /*Count num nodes for each rank (only used for preallocation)*/
        Index num_volume_elements_loc{0};
        for (Index i : part)
            if (i == i_part)
                num_volume_elements_loc++;
        volume_elements_loc.reserve(num_volume_elements_loc, MAX_NODES_VOLUME_ELEMENT);

        nodes_loc.reserve(1.5 * nodes_glob.size() / NF_MPI::get_size()); /*Assuming equal partition times a factor as estimate for number of nodes for each processor*/

        /*--------------------------------------------------------------------
        Looping over volume elements to copy from global elements and nodes
        to local elements and nodes.
        --------------------------------------------------------------------*/
        for (Index i{0}; i < volume_elements_glob.size(); i++)
        {
            if (part[i] == i_part)
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
                    if (nID_glob_to_loc.count(node_id_glob) == 0)
                    {
                        nodes_loc.emplace_back(nodes_glob[node_id_glob]);
                        Index node_id_loc = nodes_loc.size();
                        nID_glob_to_loc.emplace(node_id_glob, node_id_loc);
                        nID_loc_to_glob.emplace(node_id_loc, node_id_glob);
                    }
                }
                /*--------------------------------------------------------------------
                Copies an element from the global to the local domain, and renumbers
                the node indices from global to local values.
                --------------------------------------------------------------------*/

                assert(eID_glob_to_loc.count(i) == 0);
                eID_glob_to_loc.emplace(volume_elements_loc.size());

                volume_elements_loc.add_element_local(e_type, element, nID_glob_to_loc);
            }
        }

        /*--------------------------------------------------------------------
        Looping over patches and copying the patch elements and the patch name
        from global to local.
        --------------------------------------------------------------------*/

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
                    Index nID_glob = element[k];
                    if (nID_glob_to_loc.count(nID_glob) != 1)
                    {
                        shared_face = false;
                        break;
                    }
                }
                if (shared_face)
                    boundary_elements_loc.add_element_local(e_type, element, nID_glob_to_loc);
            }
            boundary_elements_loc.shrink_to_fit();
        }
        element_patches_loc.shrink_to_fit();
        volume_elements_loc.shrink_to_fit();
        nodes_loc.shrink_to_fit();
    }

    void GridCreator::create_FV_grid_local(const map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                           const GridUtils &grid_utils,
                                           Index rank_loc,
                                           unique_ptr<FV_Grid> &FV_Grid_loc,
                                           PrimalGrid &primal_grid_loc,
                                           map<FaceElement, InternalGhostData> const &internal_boundary_faces_glob)
    {
        const Vector<Vec3> &nodes_loc = primal_grid_loc.get_nodes();
        const Elements &volume_elements_loc = primal_grid_loc.get_volume_elements();
        Elements &face_elements_loc = primal_grid_loc.get_face_elements();
        assert(face_elements_loc.size() == 0);
        const Vector<ElementPatch> &element_patches_loc = primal_grid_loc.get_element_patches();

        Cells cells_loc;
        Faces faces_loc;
        faces_loc.reserve(1.5 * faces_to_cells_glob.size() / NF_MPI::get_size()); // Guesstimate
        Vector<Patch> patches_loc;

        /*--------------------------------------------------------------------
        Step 1: Create all cells from elements.
        --------------------------------------------------------------------*/

        cells_loc.resize(volume_elements_loc.size() + FV_Grid::find_N_GHOST_cells(element_patches_loc));

        /*--------------------------------------------------------------------
        Adding internal faces first
        --------------------------------------------------------------------*/
        for (const auto &p : faces_to_cells_glob)
        {
            Index u_glob = p.second.first;
            Index v_glob = p.second.second;
            const FaceElement &fe_glob = p.first;

            bool u_is_in_domain = grid_utils.get_partID_from_global_eID(u_glob) == i_part;
            bool v_is_in_domain = grid_utils.get_partID_from_global_eID(v_glob) == i_part;

            if (u_is_in_domain && v_is_in_domain)
            {
                Index u_loc = grid_utils.get_local_element_index(u_glob);
                Index v_loc = grid_utils.get_local_element_index(v_glob);

                Index i_loc = min(u_loc, v_loc);
                Index j_loc = max(u_loc, v_loc);

                assert(i_loc != j_loc && i_loc != NOT_ASSIGNED && j_loc != NOT_ASSIGNED);

                faces_loc.cell_indices.emplace_back(i_loc, j_loc);

                face_elements_loc.add_element_local(fe_glob.e_type,
                                                    fe_glob.nodes.data(),
                                                    grid_utils.get_map_nodeID_glob_to_loc());
            }
        }

        /*--------------------------------------------------------------------
        Step 4: Create the faces at the boundaries and ghost cells.
        Also add the face elements at the boundaries.
        --------------------------------------------------------------------*/
        cout << "Creating ghost cells..\n";
        Index j_ghost = volume_elements_loc.size();
        for (const auto &element_patch : element_patches_loc)
        {
            Patch p;
            p.boundary_type = config.get_boundary_type(element_patch.patch_name);
            p.FIRST_FACE = faces_loc.size();
            p.N_FACES = element_patch.boundary_elements.size();
            patches_loc.push_back(p);

            const Elements &surface_elements_loc = element_patch.boundary_elements;
            for (Index ij{0}; ij < surface_elements_loc.size(); ij++)
            {
                ElementType e_type = surface_elements_loc.get_element_type(ij);
                FaceElement face_element{e_type, surface_elements_loc.get_element_nodes(ij)};
                face_element.loc_to_glob(grid_utils.get_map_nodeID_loc_to_glob()); // Switching from local to global indices
                assert(faces_to_cells_glob.at(face_element).second == NOT_ASSIGNED);
                Index i_domain = faces_to_cells_glob.at(face_element).first;
                assert(i_domain < j_ghost);
                faces_loc.cell_indices.emplace_back(i_domain, j_ghost);
                face_elements_loc.add_element(e_type, face_element.nodes.data());
                j_ghost++;
            }
        }
        assert(face_elements_loc.size() == faces_loc.size());

        /*--------------------------------------------------------------------
        Step 5: Assign interior ghost cells between domains.
        Now all faces except the ones arising from interprocessor domain
        boundaries should be accounted for. Adding each domain sequentially so
        that faces are grouped together.
        --------------------------------------------------------------------*/
        for (Index rank_neigbour{0}; rank_neigbour < NF_MPI::get_size(); rank_neigbour++)
        {
            for (const auto &p : internal_boundary_faces_glob)
            {

                const InternalGhostData &gd = p.second;

                if ((gd.rank_a == rank_loc && rank_neigbour == gd.rank_b) ||
                    (gd.rank_b == rank_loc && rank_neigbour == gd.rank_a))
                {
                    Index i_loc, j_loc;
                    assert(rank_loc == gd.rank_a || rank_loc == gd.rank_b);
                    assert(rank_neigbour == gd.rank_a || rank_neigbour == gd.rank_b);
                    if (rank_loc == gd.rank_a)
                    {
                        i_loc = grid_utils.get_local_element_index(gd.cID_a);
                        j_loc = grid_utils.get_local_element_index(gd.cID_b);
                    }
                    else
                    {
                        i_loc = grid_utils.get_local_element_index(gd.cID_b);
                        j_loc = grid_utils.get_local_element_index(gd.cID_a);
                    }
                    faces_loc.cell_indices.emplace_back(i_loc, j_loc);
                    const FaceElement &fe_glob = p.first;
                    face_elements_loc.add_element_local(fe_glob.e_type, fe_glob.nodes.data(),
                                                        grid_utils.get_map_nodeID_glob_to_loc());
                }
            }
        }
        /*-------------------------------------------------------------------
            Set some grid metrics in the config object
        --------------------------------------------------------------------*/
        Index N_NODES = nodes.size();
        Index N_INTERIOR_CELLS = volume_elements.size();
        Index N_TOTAL_CELLS = cells.size();
        Index N_INTERIOR_FACES = faces.size() - find_N_GHOST_cells(element_patches);
        Index N_TOTAL_FACES = faces.size();
        config.set_grid_metrics(N_NODES, N_INTERIOR_CELLS, N_TOTAL_CELLS, N_INTERIOR_FACES, N_TOTAL_FACES);

        /*--------------------------------------------------------------------
        Sort faces so that the interior faces appear first with the owner index
        i allways being less than neigbour index j. The same logic is applied
        patch-wise to the boundaries
        --------------------------------------------------------------------*/

        /*All Vectors within faces and cells must have the correct size before reordering*/
        faces.resize_geometry_properties();

        reorder_faces(config, face_elements);

        /*--------------------------------------------------------------------
        Ensuring that no unneccessary memory isn't used
        --------------------------------------------------------------------*/
        primal_grid.element_patches.clear(); // element patches no longer needed.
        assert(primal_grid.element_patches.empty());
        assert(faces.cell_indices.capacity() == N_TOTAL_FACES);
        assert(faces.face_normals.capacity() == N_TOTAL_FACES);
        assert(faces.centroid_to_face_i.capacity() == N_TOTAL_FACES);
        assert(faces.centroid_to_face_j.capacity() == N_TOTAL_FACES);
        assert(cells.volumes.capacity() == N_TOTAL_CELLS);
        assert(cells.centroids.capacity() == N_TOTAL_CELLS);

        cout << "Computational grid has been created.\n";
    }

    void GridCreator::create_partitioned_grids(Config &config,
                                               unique_ptr<PrimalGrid> &primal_grid,
                                               unique_ptr<FV_Grid> &FV_grid)
    {

        const ShortIndex num_procs = NF_MPI::get_size();
        const ShortIndex rank = NF_MPI::get_rank();

        if (rank == 0)
        {
            auto primal_grid_glob = make_unique<PrimalGrid>(config);

            map<FaceElement, pair<Index, long int>> faces_to_cells_glob;
            Elements face_elements_glob;

            const Vector<Index> part = NF_METIS::calc_element_partition(*primal_grid_glob, num_procs);

            map<FaceElement, InternalGhostData> internal_ghost_faces;

            construct_global_face_entities(primal_grid_glob->get_volume_elements(),
                                           primal_grid_glob->get_element_patches(),
                                           faces_to_cells_glob,
                                           face_elements_glob,
                                           internal_ghost_faces,
                                           part);

            map<Index, Index> eID_glob_to_loc;
            Vector<map<Index, Index>> nID_glob_to_loc_vec(num_procs);
            Vector<map<Index, Index>> nID_loc_to_glob_vec(num_procs);
            Vector<PrimalGrid> primal_grids_loc;

            for (Index i_part{0}; i_part < num_procs; i_part++)
            {
                Elements volume_elements_loc;
                Vector<Vec3> nodes_loc;
                Vector<ElementPatch> element_patches_loc;

                create_local_volume_entities(primal_grid_glob->get_volume_elements(),
                                             primal_grid_glob->get_nodes(),
                                             primal_grid_glob->get_element_patches(),
                                             i_part,
                                             part,
                                             volume_elements_loc,
                                             nodes_loc,
                                             element_patches_loc,
                                             nID_glob_to_loc_vec[i_part],
                                             nID_loc_to_glob_vec[i_part],
                                             eID_glob_to_loc);

                primal_grids_loc.emplace_back(nodes_loc, volume_elements_loc, element_patches_loc);
            }

            for (Index i_part{0}; i_part < num_procs; i_part++)
            {
                GridUtils utils{part,
                                eID_glob_to_loc,
                                nID_glob_to_loc_vec[i_part], nID_loc_to_glob_vec[i_part]};
                unique_ptr<FV_Grid> FV_grid_loc;
                create_FV_grid_local(faces_to_cells_glob,
                                     utils,
                                     i_part,
                                     primal_grids_loc[i_part],
                                     internal_ghost_faces);
            }
        }
    }

    //     void create_partitioned_grids(Config &config,
    //                                   unique_ptr<geometry::PrimalGrid> &primal_grid,
    //                                   unique_ptr<geometry::FV_Grid> &FV_grid)
    //     {
    //         using namespace geometry;
    //         const ShortIndex num_procs = NF_MPI::get_size();
    //         const ShortIndex rank = NF_MPI::get_rank();

    //         /*--------------------------------------------------------------------
    //         Rank 0 creates a global grid.
    //         --------------------------------------------------------------------*/

    //         if (rank == 0)
    //         {
    //             auto primal_grid_glob = make_unique<PrimalGrid>(config);
    //             auto FV_Grid_glob = make_unique<FV_Grid>(config, primal_grid_glob);

    //             Vector<Index> volume_element_partition = NF_METIS::calc_element_partition(*primal_grid_glob, num_procs);

    //             // Given that I have:
    //             map<FaceElement, pair<Index, Index>> face_to_glob_cells;

    //             Vector<map<FaceElement, pair<Index, Index>>> face_to_loc_cells_vec(num_procs);

    //             for (const auto &kv : face_to_glob_cells)
    //             {
    //                 const FaceElement &face = kv.first;
    //                 Index i_glob = kv.second.first;
    //                 Index j_glob = kv.second.second;
    //                 ShortIndex proc_i = volume_element_partition[i_glob];
    //                 ShortIndex proc_j = volume_element_partition[j_glob];

    //                 if (face_to_loc_cells_vec[proc_i].count(face) == 0)
    //                     face_to_loc_cells_vec[proc_i].emplace(face, {})
    //             }

    //             /*--------------------------------------------------------------------
    //             If there is only one processor, assign the global grid as the only
    //             primal grid, construct finite volume grid, and return.
    //             --------------------------------------------------------------------*/
    //             if (num_procs == 1)
    //             {
    //                 primal_grid = move(primal_grid_glob);
    //                 FV_grid = make_unique<geometry::FV_Grid>(config, *primal_grid);
    //                 return;
    //             }

    //             /*--------------------------------------------------------------------
    //             If multiple processors are used, the domain must be partitioned into several
    //             pieces. This is handled by METIS.
    //             The return value 'volume_element_partition', is a vector of indices
    //             from 0 to n_procs-1 marking which processor each volume element
    //             belongs to.
    //             --------------------------------------------------------------------*/
    //             Vector<Index> volume_element_partition = NF_METIS::calc_element_partition(*primal_grid_glob, num_procs);

    //             const Elements &volume_elements_glob = primal_grid_glob->get_volume_elements();
    //             const Vector<Vec3> &nodes_glob = primal_grid_glob->get_nodes();
    //             assert(volume_element_partition.size() == volume_elements_glob.size());

    //             /*Count num nodes for each rank (only used for preallocation)*/
    //             Index num_volume_elements_loc{0};
    //             for (Index i : volume_element_partition)
    //                 if (i == rank)
    //                     num_volume_elements_loc++;

    //             Elements volume_elements_loc;
    //             volume_elements_loc.reserve(num_volume_elements_loc, MAX_NODES_VOLUME_ELEMENT);
    //             Vector<Vec3> nodes_loc;
    //             nodes_loc.reserve(1.5 * nodes_glob.size() / num_procs); /*Assuming equal partition times a factor as estimate for number of nodes for each processor*/

    //             map<Index, Index> glob_to_loc_node_id;

    //             /*--------------------------------------------------------------------
    //             Looping over volume elements to copy from global elements and nodes
    //             to local elements and nodes.
    //             --------------------------------------------------------------------*/
    //             for (Index i{0}; i < volume_elements_glob.size(); i++)
    //             {
    //                 if (volume_element_partition[i] == rank)
    //                 {
    //                     const Index *element = volume_elements_glob.get_element_nodes(i);
    //                     ShortIndex num_nodes = volume_elements_glob.get_n_element_nodes(i);
    //                     ElementType e_type = volume_elements_glob.get_element_type(i);

    //                     /*--------------------------------------------------------------------
    //                     Loops over all global node indices of an element. If a new node index is
    //                     discovered, the point corresponding to the node index is added and the
    //                     a mapping from global indices are added.
    //                     --------------------------------------------------------------------*/
    //                     for (ShortIndex k{0}; k < num_nodes; k++)
    //                     {
    //                         Index node_id_glob = element[k];
    //                         if (glob_to_loc_node_id.count(node_id_glob) == 0)
    //                         {
    //                             nodes_loc.emplace_back(nodes_glob[node_id_glob]);
    //                             Index node_id_loc = nodes_loc.size();
    //                             glob_to_loc_node_id.emplace(node_id_glob, node_id_loc);
    //                         }
    //                     }
    //                     /*--------------------------------------------------------------------
    //                     Copies an element from the global to the local domain, and renumbers
    //                     the node indices from global to local values.
    //                     --------------------------------------------------------------------*/
    //                     volume_elements_loc.add_element_local(e_type, element, glob_to_loc_node_id);
    //                 }
    //             }

    //             /*--------------------------------------------------------------------
    //             Looping over patches and copying the patch elements and the patch name
    //             from global to local.
    //             --------------------------------------------------------------------*/
    //             Vector<ElementPatch> element_patches_loc;

    //             const Vector<ElementPatch> &element_patches_glob = primal_grid->get_element_patches();

    //             for (const auto &element_patch_glob : element_patches_glob)
    //             {
    //                 const Elements &boundary_elements_glob = element_patch_glob.boundary_elements;
    //                 const string &patch_name_glob = element_patch_glob.patch_name;
    //                 element_patches_loc.emplace_back();
    //                 element_patches_loc.back().patch_name = patch_name_glob;
    //                 Elements &boundary_elements_loc = element_patches_loc.back().boundary_elements;
    //                 boundary_elements_loc.reserve(boundary_elements_glob.size(), MAX_NODES_FACE_ELEMENT);

    //                 /*--------------------------------------------------------------------
    //                 Looping over all elements of a patch
    //                 --------------------------------------------------------------------*/
    //                 for (Index i{0}; i < boundary_elements_glob.size(); i++)
    //                 {
    //                     ElementType e_type = boundary_elements_glob.get_element_type(i);
    //                     const Index *element = boundary_elements_glob.get_element_nodes(i);
    //                     ShortIndex num_nodes = boundary_elements_glob.get_n_element_nodes(i);

    //                     /*--------------------------------------------------------------------
    //                     Using the global to local node index map (from the volume element loop)
    //                     to check wether a boundary face belongs to the local domain. Keep in mind
    //                     that this map only contains the nodes belonging to the local domain.
    //                     If all global nodes are common, the global face has to be part of the
    //                     local domain. In that case, the boundary element is added. Interfaces
    //                     between local domains are not handled in this step.
    //                     --------------------------------------------------------------------*/
    //                     bool shared_face = true;
    //                     for (ShortIndex k{0}; k < num_nodes; k++)
    //                     {
    //                         Index node_id_glob = element[k];
    //                         if (glob_to_loc_node_id.count(node_id_glob) != 1)
    //                         {
    //                             shared_face = false;
    //                             break;
    //                         }
    //                     }
    //                     if (shared_face)
    //                         boundary_elements_loc.add_element_local(e_type, element, glob_to_loc_node_id);
    //                 }
    //                 boundary_elements_loc.shrink_to_fit();
    //             }
    //             element_patches_loc.shrink_to_fit();
    //             volume_elements_loc.shrink_to_fit();
    //             nodes_loc.shrink_to_fit();

    //             primal_grid = make_unique<PrimalGrid>(nodes_loc, volume_elements_loc, element_patches_loc);
    //         }
    //     }
    // }
    namespace NF_METIS
    {
        Vector<Index> calc_element_partition(PrimalGrid &primal_grid,
                                             Index n_partitions)
        {
            static_assert(sizeof(idx_t) == sizeof(Index));

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