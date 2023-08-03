#include "../../include/parallelization/GridCreator.hpp"

/*--------------------------------------------------------------------
This function reads the mesh file and creates both a primal grid
and a finite volume grid. If multiple processors are used, Metis
is used to partition the mesh.
--------------------------------------------------------------------*/

namespace geometry
{

    void GridCreator::create_partitioned_grids(Config &config,
                                               unique_ptr<PrimalGrid> &primal_grid_glob,
                                               unique_ptr<PrimalGrid> &primal_grid,
                                               unique_ptr<FV_Grid> &FV_grid)
    {

        const ShortIndex num_procs = NF_MPI::get_size();
        const ShortIndex rank = NF_MPI::get_rank();

        Vector<unique_ptr<PrimalGrid>> primal_grids_loc(num_procs);
        Vector<unique_ptr<FV_Grid>> FV_grids_loc(num_procs);

        if (rank == 0)
        {
            primal_grid_glob = make_unique<PrimalGrid>(config);
            const Vector<Index> part = NF_METIS::calc_element_partition(*primal_grid_glob, num_procs);
            Vector<pair<Index, Index>> part_to_element_range;
            Vector<Index> eID_glob_to_loc;
            reorder_global_grid(part, *primal_grid_glob, part_to_element_range, eID_glob_to_loc);

            map<FaceElement, pair<Index, long int>> faces_to_cells_glob;
            Elements face_elements_glob;

            map<FaceElement, GhostDataPartition> internal_ghost_faces;

            create_global_face_entities(primal_grid_glob->get_vol_elements(),
                                        primal_grid_glob->get_element_patches(),
                                        faces_to_cells_glob,
                                        face_elements_glob,
                                        internal_ghost_faces,
                                        part);

            Vector<map<Index, Index>> nID_glob_to_loc_vec(num_procs);
            Vector<map<Index, Index>> nID_loc_to_glob_vec(num_procs);

            /*--------------------------------------------------------------------
            Build local primal grids from global primal grid (and some addressing
            structures)
            --------------------------------------------------------------------*/

            for (Index r_loc{0}; r_loc < num_procs; r_loc++)
            {
                Elements vol_elements_loc;
                Vector<Vec3> nodes_loc;
                Vector<ElementPatch> element_patches_loc;

                create_primal_grid_local(primal_grid_glob->get_vol_elements(),
                                         primal_grid_glob->get_nodes(),
                                         primal_grid_glob->get_element_patches(),
                                         r_loc,
                                         part_to_element_range,
                                         part,
                                         primal_grids_loc[r_loc],
                                         eID_glob_to_loc,
                                         nID_glob_to_loc_vec[r_loc],
                                         nID_loc_to_glob_vec[r_loc]);
            }
            /*--------------------------------------------------------------------
            Build local FV_Grids from local primal grids
            --------------------------------------------------------------------*/

            for (Index r_loc{0}; r_loc < num_procs; r_loc++)
            {
                PartitionUtils utils{part,
                                     part_to_element_range,
                                     eID_glob_to_loc,
                                     nID_glob_to_loc_vec[r_loc], nID_loc_to_glob_vec[r_loc]};
                create_FV_grid_local(config,
                                     faces_to_cells_glob,
                                     utils,
                                     r_loc,
                                     FV_grids_loc[r_loc],
                                     *primal_grids_loc[r_loc],
                                     internal_ghost_faces);
            }
            assert(primal_grids_loc.size() == num_procs && FV_grids_loc.size() == num_procs);

            /*Setting global grid metrics in config for rank 0*/
            Index N_NODES_GLOB = primal_grid_glob->get_nodes().size();
            Index N_INTERIOR_CELLS_GLOB = primal_grid_glob->get_vol_elements().size();
            Index n_ghost_tot = primal_grid_glob->find_num_ghost_external();
            Index N_TOTAL_CELLS_GLOB = primal_grid_glob->get_vol_elements().size() + n_ghost_tot;
            Index N_TOTAL_FACES_GLOB = primal_grid_glob->get_face_elements().size();
            Index N_INTERIOR_FACES_GLOB = N_TOTAL_FACES_GLOB - n_ghost_tot;
            config.set_grid_metrics_global(N_NODES_GLOB, N_INTERIOR_CELLS_GLOB,
                                           N_TOTAL_CELLS_GLOB, N_INTERIOR_FACES_GLOB, N_TOTAL_FACES_GLOB);
        }
        send_recv_grids(config, primal_grids_loc, FV_grids_loc, primal_grid, FV_grid);
        set_config_grid_data_local(config, primal_grid, FV_grid);
    }

    /*Reordering volume elements of primal grid so that the elements within each partition are clustered
    together*/
    void GridCreator::reorder_global_grid(const Vector<Index> &part,
                                          PrimalGrid &primal_grid,
                                          Vector<pair<Index, Index>> &part_to_element_range,
                                          Vector<Index> &eID_glob_to_loc)
    {
        assert(NF_MPI::get_rank() == 0);
        ShortIndex num_procs = NF_MPI::get_size();
        part_to_element_range.resize(num_procs);

        Elements &vol_elements_old = primal_grid.get_vol_elements();
        Index n_elem = vol_elements_old.size();
        eID_glob_to_loc.resize(n_elem);
        Elements vol_elements_new;
        vol_elements_new.reserve(n_elem, MAX_NODES_VOLUME_ELEMENT);
        for (ShortIndex r_loc{0}; r_loc < num_procs; r_loc++)
        {
            /*This loop can be made faster (outermost loop removed) by first creating part_to_element_range
        (count occurences of each rank) and then using the intervals to insert elements in vol_elements_new*/

            bool first_element_found = false;
            for (Index i{0}; i < n_elem; i++)
            {
                if (part[i] == r_loc)
                {
                    if (!first_element_found)
                        part_to_element_range[r_loc].first = i;
                    first_element_found = true;
                    part_to_element_range[r_loc].second = i + 1;
                    vol_elements_new.add_element(vol_elements_old.get_element_type(i), vol_elements_old.get_element_nodes(i));
                    eID_glob_to_loc[i] = //???FIX THIS
                }
            }
        }
        assert(vol_elements_new.size() == vol_elements_old.size());
        assert(eID_glob_to_loc.size() == vol_elements_old.size());
        vol_elements_new.shrink_to_fit();
        vol_elements_old = vol_elements_new;
    }

    void GridCreator::create_global_face_entities(const Elements &vol_elements_glob,
                                                  const Vector<ElementPatch> &element_patches,
                                                  map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                                  Elements &face_elements_glob,
                                                  map<FaceElement, GhostDataPartition> &internal_boundary_faces,
                                                  const Vector<Index> &part)
    {
        assert(NF_MPI::get_rank() == 0);

        assert(faces_to_cells_glob.size() == 0);

        for (Index cell_index{0}; cell_index < vol_elements_glob.size(); cell_index++)
        {
            ElementType volume_element_type = vol_elements_glob.get_element_type(cell_index);
            assert(is_volume_element(volume_element_type));
            for (ShortIndex k{0}; k < get_num_faces_volume_element(volume_element_type); k++)
            {
                FaceElement face_element = get_face_element_k_of_volume_element(volume_element_type,
                                                                                vol_elements_glob.get_element_nodes(cell_index),
                                                                                k);

                assert(faces_to_cells_glob.count(face_element) <= 1);
                if (faces_to_cells_glob.count(face_element) == 0)
                {
                    // Discovered a new face
                    faces_to_cells_glob.emplace(face_element, pair{cell_index, CELL_NOT_ASSIGNED});
                }
                else
                {
                    assert(faces_to_cells_glob.at(face_element).second == CELL_NOT_ASSIGNED);
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
            if (cell_j != CELL_NOT_ASSIGNED)
            { // Only add internal faces
                assert(cell_i < cell_j);
                face_elements_glob.add_element(face.first.e_type, face.first.nodes.data());

                /*Marking partition boundaries*/
                if (part[cell_i] != part[cell_j])
                {
                    assert(internal_boundary_faces.count(face.first) == 0);

                    internal_boundary_faces.emplace(face.first, GhostDataPartition{cell_i,
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

        // Index j_ghost = vol_elements_glob.size();
        for (const auto &element_PatchExt : element_patches)
        {

            const Elements &surface_elements = element_PatchExt.boundary_elements;
            for (Index ij{0}; ij < surface_elements.size(); ij++)
            {
                ElementType e_type = surface_elements.get_element_type(ij);
                FaceElement face_element{e_type, surface_elements.get_element_nodes(ij)};
                assert(faces_to_cells_glob.at(face_element).second == CELL_NOT_ASSIGNED);
                // faces_to_cells.at(face_element).second = j_ghost;
                Index i_domain = faces_to_cells_glob.at(face_element).first;
                face_elements_glob.add_element(e_type, face_element.nodes.data());
                // j_ghost++;
            }
        }
        assert(face_elements_glob.size() == faces_to_cells_glob.size());
    }

    void GridCreator::create_primal_grid_local(const Elements &vol_elements_glob,
                                               const Vector<Vec3> &nodes_glob,
                                               const Vector<ElementPatch> &element_patches_glob,
                                               Index r_loc,
                                               const Vector<pair<Index, Index>> &part_to_element_range,
                                               const Vector<Index> &part,
                                               unique_ptr<PrimalGrid> &primal_grid_loc,
                                               Vector<Index> &eID_glob_to_loc,
                                               map<Index, Index> &nID_glob_to_loc,
                                               map<Index, Index> &nID_loc_to_glob)
    {

        assert(NF_MPI::get_rank() == 0);

        Vector<Vec3> nodes_loc;
        Elements vol_elements_loc;
        Vector<ElementPatch> element_patches_loc;

        assert(part.size() == vol_elements_glob.size());
        assert(vol_elements_loc.size() == 0);
        assert(nodes_loc.size() == 0);
        assert(element_patches_loc.size() == 0);
        assert(nID_glob_to_loc.size() == 0);

        Index num_vol_elements_loc = part_to_element_range[r_loc].second - part_to_element_range[r_loc].first;
#ifndef NDEBUG
        Index num_vol_elements_loc_test{0};
        for (Index i : part)
            if (i == r_loc)
                num_vol_elements_loc_test++;
        assert(num_vol_elements_loc == num_vol_elements_loc_test);
#endif

        vol_elements_loc.reserve(num_vol_elements_loc, MAX_NODES_VOLUME_ELEMENT);

        nodes_loc.reserve(1.5 * nodes_glob.size() / NF_MPI::get_size()); /*Assuming equal partition times a factor as estimate for number of nodes for each processor*/

        /*--------------------------------------------------------------------
        Looping over volume elements to copy from global elements and nodes
        to local elements and nodes.
        --------------------------------------------------------------------*/
        Index loc_e_begin = part_to_element_range[r_loc].first;
        Index loc_e_end = part_to_element_range[r_loc].second;

#ifndef NDEBUG
        if (r_loc == 0)
            assert(loc_e_begin == 0);
        else if (r_loc == NF_MPI::get_size() - 1)
            assert(loc_e_end == vol_elements_glob.size());
#endif

        for (Index i{loc_e_begin}; i < loc_e_end; i++)
        {
            assert(part[i] == r_loc);

            const Index *element = vol_elements_glob.get_element_nodes(i);
            ShortIndex num_nodes = vol_elements_glob.get_n_element_nodes(i);
            ElementType e_type = vol_elements_glob.get_element_type(i);

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
            vol_elements_loc.add_element_local(e_type, element, nID_glob_to_loc);
            assert(vol_elements_loc.size() == eID_glob_to_loc[i]);
        }

        /*--------------------------------------------------------------------
        Looping over patches and copying the PatchExt elements and the PatchExt name
        from global to local.
        --------------------------------------------------------------------*/

        for (const auto &element_PatchExt_glob : element_patches_glob)
        {
            const Elements &boundary_elements_glob = element_PatchExt_glob.boundary_elements;
            const string &PatchExt_name_glob = element_PatchExt_glob.PatchExt_name;
            element_patches_loc.emplace_back();
            element_patches_loc.back().PatchExt_name = PatchExt_name_glob;
            Elements &boundary_elements_loc = element_patches_loc.back().boundary_elements;
            boundary_elements_loc.reserve(boundary_elements_glob.size(), MAX_NODES_FACE_ELEMENT);

            /*--------------------------------------------------------------------
            Looping over all elements of a PatchExt
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
        vol_elements_loc.shrink_to_fit();
        nodes_loc.shrink_to_fit();
        primal_grid_loc = make_unique<PrimalGrid>(nodes_loc, vol_elements_loc, element_patches_loc, loc_e_begin);
    }

    void GridCreator::create_FV_grid_local(const Config &config,
                                           const map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                           const PartitionUtils &utils,
                                           Index r_loc,
                                           unique_ptr<FV_Grid> &FV_grid_loc,
                                           PrimalGrid &primal_grid_loc,
                                           map<FaceElement, GhostDataPartition> const &internal_boundary_faces_glob)
    {
        const Vector<Vec3> &nodes_loc = primal_grid_loc.get_nodes();
        const Elements &vol_elements_loc = primal_grid_loc.get_vol_elements();
        Elements &face_elements_loc = primal_grid_loc.get_face_elements();
        assert(face_elements_loc.size() == 0);
        const Vector<ElementPatch> &element_patches_loc = primal_grid_loc.get_element_patches();

        Cells cells_loc;
        Faces faces_loc;
        faces_loc.reserve(1.5 * faces_to_cells_glob.size() / NF_MPI::get_size()); // Guesstimate
        Vector<PatchExt> patches_loc;
        Vector<PartitionPatchExt> part_patches_loc;

        /*--------------------------------------------------------------------
        Adding internal faces first
        --------------------------------------------------------------------*/
        for (const auto &p : faces_to_cells_glob)
        {
            Index u_glob = p.second.first;
            Index v_glob = p.second.second;
            const FaceElement &fe_glob = p.first;

            bool u_is_in_domain = utils.get_partID_from_global_eID(u_glob) == r_loc;
            bool v_is_in_domain = utils.get_partID_from_global_eID(v_glob) == r_loc;

            if (u_is_in_domain && v_is_in_domain)
            {
                Index u_loc = utils.get_local_element_index(u_glob);
                Index v_loc = utils.get_local_element_index(v_glob);

                Index i_loc = min(u_loc, v_loc);
                Index j_loc = max(u_loc, v_loc);

                assert(i_loc != j_loc && i_loc != CELL_NOT_ASSIGNED && j_loc != CELL_NOT_ASSIGNED);

                faces_loc.cell_indices.emplace_back(i_loc, j_loc);

                face_elements_loc.add_element_local(fe_glob.e_type,
                                                    fe_glob.nodes.data(),
                                                    utils.get_map_nodeID_glob_to_loc());
            }
        }

        /*--------------------------------------------------------------------
        Assign interior ghost cells between domains.
        Now all faces except the ones arising from interprocessor domain
        boundaries should be accounted for. Adding each domain sequentially so
        that faces are grouped together.
        --------------------------------------------------------------------*/
        Index num_ghost_part{0};
        for (Index rank_neigbour{0}; rank_neigbour < NF_MPI::get_size(); rank_neigbour++)
        {
            Index num_PatchExt_faces = 0;
            for (const auto &p : internal_boundary_faces_glob)
            {

                const GhostDataPartition &gd = p.second;

                if ((gd.rank_a == r_loc && rank_neigbour == gd.rank_b) ||
                    (gd.rank_b == r_loc && rank_neigbour == gd.rank_a))
                {
                    num_PatchExt_faces++;
                    num_ghost_part++;
                    Index i_loc, j_loc;
                    assert(r_loc == gd.rank_a || r_loc == gd.rank_b);
                    assert(rank_neigbour == gd.rank_a || rank_neigbour == gd.rank_b);
                    if (r_loc == gd.rank_a)
                    {
                        i_loc = utils.get_local_element_index(gd.cID_a);
                        j_loc = utils.get_local_element_index(gd.cID_b);
                    }
                    else
                    {
                        i_loc = utils.get_local_element_index(gd.cID_b);
                        j_loc = utils.get_local_element_index(gd.cID_a);
                    }
                    faces_loc.cell_indices.emplace_back(i_loc, j_loc);
                    const FaceElement &fe_glob = p.first;
                    face_elements_loc.add_element_local(fe_glob.e_type, fe_glob.nodes.data(),
                                                        utils.get_map_nodeID_glob_to_loc());
                }
            }
            /*--------------------------------------------------------------------
            For now I will just create patches for all neighbour ranks for the local
            rank regardless of they share faces or not. This is to avoid potential
            communication locking further down the pipeline. I should look into this
            and change it in the future, so that only shared patches are added.
            --------------------------------------------------------------------*/
            PartitionPatchExt part_PatchExt;
            part_PatchExt.FIRST_FACE = faces_loc.size();
            part_PatchExt.N_FACES = num_PatchExt_faces;
            part_PatchExt.rank_neighbour = rank_neigbour;
            part_patches_loc.emplace_back(part_PatchExt);
        }

        /*--------------------------------------------------------------------
        Create the faces at the boundaries and ghost cells.
        Also add the face elements at the boundaries.
        --------------------------------------------------------------------*/

        Index j_ghost = vol_elements_loc.size();
        for (const auto &element_PatchExt : element_patches_loc)
        {
            PatchExt p;
            p.boundary_type = config.get_boundary_type(element_PatchExt.PatchExt_name);
            p.FIRST_FACE = faces_loc.size();
            p.N_FACES = element_PatchExt.boundary_elements.size();
            patches_loc.push_back(p);

            const Elements &surface_elements_loc = element_PatchExt.boundary_elements;
            for (Index ij{0}; ij < surface_elements_loc.size(); ij++)
            {
                ElementType e_type = surface_elements_loc.get_element_type(ij);
                FaceElement face_element{e_type, surface_elements_loc.get_element_nodes(ij)};
                face_element.loc_to_glob(utils.get_map_nodeID_loc_to_glob()); // Switching from local to global indices
                assert(faces_to_cells_glob.at(face_element).second == CELL_NOT_ASSIGNED);
                Index i_domain = faces_to_cells_glob.at(face_element).first;
                assert(i_domain < j_ghost);
                faces_loc.cell_indices.emplace_back(i_domain, j_ghost);
                face_elements_loc.add_element(e_type, face_element.nodes.data());
                j_ghost++;
            }
        }
        assert(face_elements_loc.size() == faces_loc.size());

        Index tot_ghost_loc = primal_grid_loc.find_num_ghost_external() + num_ghost_part;

        /*--------------------------------------------------------------------
        Create all cells.
        --------------------------------------------------------------------*/
        cells_loc.resize(vol_elements_loc.size() + tot_ghost_loc);

        /*-------------------------------------------------------------------
            Set some grid metrics in the config object
        --------------------------------------------------------------------*/

        faces_loc.resize_geometry_properties();

        Index num_interior_faces_loc = faces_loc.size() - tot_ghost_loc;
        reorder_face_enitities(num_interior_faces_loc, part_patches_loc, patches_loc, faces_loc, face_elements_loc);

        /*--------------------------------------------------------------------
        Ensuring that no unneccessary memory isn't used
        --------------------------------------------------------------------*/
        primal_grid_loc.element_patches.clear(); // element patches no longer needed.
        assert(faces_loc.cell_indices.capacity() == faces_loc.size());
        assert(faces_loc.face_normals.capacity() == faces_loc.size());
        assert(faces_loc.centroid_to_face_i.capacity() == faces_loc.size());
        assert(faces_loc.centroid_to_face_j.capacity() == faces_loc.size());
        assert(cells_loc.volumes.capacity() == cells_loc.size());
        assert(cells_loc.centroids.capacity() == cells_loc.size());
        assert(part_patches_loc.size() >= r_loc && part_patches_loc.size() <= NF_MPI::get_size());

        FV_grid_loc = make_unique<FV_Grid>(cells_loc, faces_loc, patches_loc, part_patches_loc);

        // cout << "Computational grid has been created.\n";
    }

    void GridCreator::reorder_face_enitities(Index num_interior_faces,
                                             const Vector<PartitionPatchExt> &partition_patches,
                                             const Vector<PatchExt> &patches,
                                             Faces &faces,
                                             Elements &face_elements)
    {
        /*Sorts the faces based on the indices of the neighbour cells. We want them to be sorted in increasing order,
        with the owner cell i being prioritized first and the neighbour cell j second. As an example,
        the ordering {(0,1), (0,3), (0,2), (1,1)} should be changed to {(0,1), (0,2), (0,3), (1,1)}
        The way the grid is constructed this is allready achieved for the owner cell i. However j might need sorting.
        The comparison operator < for Face is defined for this purpose.
        The interior and each boundary PatchExt is sorted separately. The face elements are sorted accordingly*/

        Elements face_elements_to_sort;

        faces.sort_face_entities(0, num_interior_faces, face_elements, face_elements_to_sort);

        for (const PartitionPatchExt &p_PatchExt : partition_patches)
        {
            faces.sort_face_entities(p_PatchExt.FIRST_FACE, p_PatchExt.FIRST_FACE + p_PatchExt.N_FACES, face_elements, face_elements_to_sort);
        }

        for (const PatchExt &PatchExt : patches)
        {
            faces.sort_face_entities(PatchExt.FIRST_FACE, PatchExt.FIRST_FACE + PatchExt.N_FACES, face_elements, face_elements_to_sort);
        }
        face_elements = face_elements_to_sort;
    }

    /*Uses MPI to copy all the local PrimalGrids and FV_Grids from rank 0 to the other ranks, it
    also sets various data in config objects*/
    void GridCreator::send_recv_grids(Config &config,
                                      Vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                      Vector<unique_ptr<FV_Grid>> &FV_grids_loc,
                                      unique_ptr<PrimalGrid> &primal_grid,
                                      unique_ptr<FV_Grid> &FV_grid)
    {
        ShortIndex rank = NF_MPI::get_rank();
        ShortIndex num_procs = NF_MPI::get_size();

        if (rank == 0)
        {
            primal_grid = move(primal_grids_loc[0]);
            FV_grid = move(FV_grids_loc[0]);
            for (ShortIndex r_loc{1}; r_loc < num_procs; r_loc++)
            {
                /*Send primal grid*/
                string bytes;
                serialization::serialize(bytes, *primal_grids_loc[r_loc]);
                Index num_bytes = bytes.size();
                NF_MPI::Send(&num_bytes, 1, r_loc);
                NF_MPI::Send(bytes.data(), num_bytes, r_loc);

                /*Send FV grid */
                bytes.clear();
                serialization::serialize(bytes, *FV_grids_loc[r_loc]);
                num_bytes = bytes.size();
                NF_MPI::Send(&num_bytes, 1, r_loc);
                NF_MPI::Send(bytes.data(), num_bytes, r_loc);
            }
        }
        else
        {
            /*Receive primal grid*/
            Index num_bytes;
            NF_MPI::Recv(&num_bytes, 1, 0);
            string bytes;
            bytes.resize(num_bytes);
            NF_MPI::Recv(bytes.data(), num_bytes, 0);
            serialization::deserialize(bytes, *primal_grid);

            /*Receive FV grid*/
            NF_MPI::Recv(&num_bytes, 1, 0);
            bytes.resize(num_bytes);
            NF_MPI::Recv(bytes.data(), num_bytes, 0);
            serialization::deserialize(bytes, *FV_grid);
        }

        /*Send Config global grid data*/
        Index N_NODES_GLOB, N_INTERIOR_CELLS_GLOB, N_TOTAL_CELLS_GLOB, N_TOTAL_FACES_GLOB, N_INTERIOR_FACES_GLOB;
        if (rank == 0)
        {
            N_NODES_GLOB = config.get_N_NODES_GLOB();
            N_INTERIOR_CELLS_GLOB = config.get_N_INTERIOR_CELLS_GLOB();
            N_TOTAL_CELLS_GLOB = config.get_N_TOTAL_CELLS_GLOB();
            N_INTERIOR_FACES_GLOB = config.get_N_INTERIOR_FACES_GLOB();
            N_TOTAL_FACES_GLOB = config.get_N_TOTAL_FACES_GLOB();
        }
        NF_MPI::Bcast(&N_NODES_GLOB, 1, 0);
        NF_MPI::Bcast(&N_INTERIOR_CELLS_GLOB, 1, 0);
        NF_MPI::Bcast(&N_TOTAL_CELLS_GLOB, 1, 0);
        NF_MPI::Bcast(&N_INTERIOR_FACES_GLOB, 1, 0);
        NF_MPI::Bcast(&N_TOTAL_FACES_GLOB, 1, 0);
        if (rank != 0)
            config.set_grid_metrics_global(N_NODES_GLOB,
                                           N_INTERIOR_CELLS_GLOB,
                                           N_TOTAL_CELLS_GLOB,
                                           N_INTERIOR_FACES_GLOB,
                                           N_TOTAL_FACES_GLOB);

        NF_MPI::Barrier();

        set_config_grid_data_local(config, primal_grid, FV_grid);
    }

    void GridCreator::set_config_grid_data_local(Config &config,
                                                 unique_ptr<PrimalGrid> &primal_grid,
                                                 unique_ptr<FV_Grid> &FV_grid)
    {
        Index N_NODES_LOC = primal_grid->get_nodes().size();
        Index N_INTERIOR_CELLS_LOC = primal_grid->get_vol_elements().size();
        Index N_TOTAL_CELLS_LOC = FV_grid->get_cells().size();
        Index N_INTERIOR_FACES_LOC = FV_grid->get_faces().size() - FV_grid->find_num_ghost_tot();
        Index N_TOTAL_FACES_LOC = FV_grid->get_faces().size();
        config.set_grid_metrics_local(N_NODES_LOC, N_INTERIOR_CELLS_LOC, N_TOTAL_CELLS_LOC, N_INTERIOR_FACES_LOC, N_TOTAL_FACES_LOC);
    }
}
