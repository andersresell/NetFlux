#include "../../include/parallelization/GridCreator.hpp"

namespace geometry
{ /*Main loop*/
    void GridCreator::create_partitioned_grids(Config &config,
                                               unique_ptr<PrimalGrid> &primal_grid_glob,
                                               unique_ptr<PrimalGrid> &primal_grid,
                                               unique_ptr<FV_Grid> &FV_grid,
                                               unique_ptr<PartitionComm> &part_comm)
    {
        num_part = NF_MPI::get_size();
        // num_part = 2; // dbg
        // assert(NF_MPI::get_size() == 1);

        const ShortIndex rank = NF_MPI::get_rank();

        primal_grid = make_unique<PrimalGrid>();
        FV_grid = make_unique<FV_Grid>();

        vector<unique_ptr<PrimalGrid>> primal_grids_loc;
        vector<unique_ptr<FV_Grid>> FV_grids_loc;

        if (rank == 0)
        {
            /*Step 1*/
            primal_grid_glob = make_unique<PrimalGrid>(config);
            /*Step 2*/
            const vector<ShortIndex> part = NF_METIS::calc_element_partition(*primal_grid_glob, num_part);
            Utils utils{part};

            /*Step 3*/
            reorder_global_grid(*primal_grid_glob, utils);

            /*Step 4*/
            FaceGraphGlob face_graph_glob{part};
            Elements face_elements_glob;
            create_global_face_graph_and_face_elements(config,
                                                       *primal_grid_glob,
                                                       face_graph_glob,
                                                       utils);

            /*Step 5*/
            vector<FaceGraphLoc> face_graphs_loc = create_local_face_graphs(face_graph_glob,
                                                                            utils);

            /*Step 6*/
            create_local_primal_grids(*primal_grid_glob,
                                      face_graphs_loc,
                                      utils,
                                      primal_grids_loc);

            /*Step 7*/
            create_local_FV_grids(face_graphs_loc,
                                  primal_grids_loc,
                                  utils,
                                  FV_grids_loc);

            set_global_config_data(config, primal_grid_glob, face_graph_glob);
        }
        /*Step 8*/
        send_recv_grids(config,
                        primal_grids_loc,
                        FV_grids_loc,
                        primal_grid,
                        FV_grid);

        /*Step 9*/
        FV_grid->resize_geometry_properties();
        part_comm = make_unique<PartitionComm>(config, FV_grid->faces, FV_grid->patches_interf);
        FV_grid->calc_geometry_properties(config, *primal_grid, *part_comm);

        cout << "Grid created\n";
    }

    /*Step 3*/
    void GridCreator::reorder_global_grid(PrimalGrid &primal_grid,
                                          Utils &utils)
    {
        assert(NF_MPI::get_rank() == 0);

        vector<pair<Index, Index>> part2e_range(num_part);

        Elements &vol_elements_old = primal_grid.get_vol_elements();
        Index n_elem = vol_elements_old.size();
        vector<Index> eIDglob2loc(n_elem);
        Elements vol_elements_new;
        vol_elements_new.reserve(n_elem, MAX_NODES_VOLUME_ELEMENT);
        for (ShortIndex r_loc{0}; r_loc < num_part; r_loc++)
        {
            /*This loop can be made faster (outermost loop removed) by first creating part_to_element_range
        (count occurences of each rank) and then using the intervals to insert elements in vol_elements_new*/

            bool first_element_found = false;
            for (Index i{0}; i < n_elem; i++)
            {
                if (utils.e2r(i) == r_loc)
                {
                    vol_elements_new.add_element(vol_elements_old.get_element_type(i), vol_elements_old.get_element_nodes(i));
                    Index eIDglob = vol_elements_new.size() - 1;

                    if (!first_element_found)
                    {
                        part2e_range[r_loc].first = eIDglob;
                        first_element_found = true;
                    }

                    eIDglob2loc[i] = eIDglob - part2e_range[r_loc].first;
                    // eIDglob2loc[i] = (r_loc == 0) ? eIDglob : eIDglob - part2e_range[r_loc].first;
                }
                part2e_range[r_loc].second = vol_elements_new.size();
            }
        }
        assert(vol_elements_new.size() == vol_elements_old.size());
        assert(eIDglob2loc.size() == vol_elements_old.size());
        utils.set_eIDglob2loc(move(eIDglob2loc));
        utils.set_part2e_range(move(part2e_range));
        vol_elements_new.shrink_to_fit();
        vol_elements_old = vol_elements_new;
    }

    /*Step 4*/
    void GridCreator::create_global_face_graph_and_face_elements(const Config &config,
                                                                 PrimalGrid &primal_grid,
                                                                 FaceGraphGlob &face_graph,
                                                                 Utils &utils)
    {
        Elements &face_elements = primal_grid.get_face_elements();
        assert(face_graph.size() == 0);
        assert(face_elements.size() == 0);

        const Elements &vol_elements = primal_grid.get_vol_elements();
        const vector<ElementPatch> &element_patches = primal_grid.get_element_patches();

        /*--------------------------------------------------------------------
        Associate the face nodes with its two neighboring cells.
        --------------------------------------------------------------------*/

        map<FaceElement, pair<Index, SignedIndex>> faces_to_cells;
        constexpr SignedIndex CELL_NOT_YET_ASSIGNED{-1};

        for (Index cell_index{0}; cell_index < vol_elements.size(); cell_index++)
        {
            ElementType volume_element_type = vol_elements.get_element_type(cell_index);
            assert(is_volume_element(volume_element_type));
            for (ShortIndex k{0}; k < get_num_faces_volume_element(volume_element_type); k++)
            {
                FaceElement face_element = get_face_element_k_of_volume_element(volume_element_type,
                                                                                vol_elements.get_element_nodes(cell_index),
                                                                                k);

                assert(faces_to_cells.count(face_element) <= 1);
                if (faces_to_cells.count(face_element) == 0)
                {
                    // Discovered a new face
                    faces_to_cells.emplace(face_element, pair{cell_index, CELL_NOT_YET_ASSIGNED});
                }
                else
                {
                    assert(faces_to_cells.at(face_element).second == CELL_NOT_YET_ASSIGNED);
                    // Face allready discovered by previous cell
                    faces_to_cells.at(face_element).second = cell_index;
                }
            }
        }
        face_elements.reserve(faces_to_cells.size(), MAX_NODES_FACE_ELEMENT);

        /*--------------------------------------------------------------------
        Adding internal faces to the map. (Boundary faces are added
        later, this is to get the correct grouping of patches). Face elements
        (used for calculating geometry properties) are created simultaneously
        as faces, to get the correct ordering.
        --------------------------------------------------------------------*/

        for (const auto &face : faces_to_cells)
        {
            Index cell_i = face.second.first;
            long int cell_j = face.second.second;
            if (cell_j != CELL_NOT_YET_ASSIGNED)
            { // Only add internal faces
                assert(cell_i < cell_j);
                face_graph.add_face(cell_i, cell_j);
                face_elements.add_element(face.first.e_type, face.first.nodes.data());
            }
        }
        assert(face_elements.size() == face_graph.size());

        /*--------------------------------------------------------------------
        Create the faces at the boundaries and ghost cells.
        Also add the face elements at the boundaries.
        --------------------------------------------------------------------*/

        Index j_ghost = vol_elements.size();
        for (const auto &e_patch : element_patches)
        {
            PatchBoundary p;
            p.boundary_type = config.get_boundary_type(e_patch.patch_name);
            p.FIRST_FACE = face_graph.size();
            p.N_FACES = e_patch.boundary_elements.size();
            face_graph.add_patch_bound(p);

            const Elements &surface_elements = e_patch.boundary_elements;
            for (Index ij{0}; ij < surface_elements.size(); ij++)
            {
                ElementType e_type = surface_elements.get_element_type(ij);
                FaceElement face_element{e_type, surface_elements.get_element_nodes(ij)};
                assert(faces_to_cells.at(face_element).second == CELL_NOT_YET_ASSIGNED);
                faces_to_cells.at(face_element).second = j_ghost;
                Index i_domain = faces_to_cells.at(face_element).first;
                assert(i_domain < j_ghost);
                // face_graph.add_face(i_domain,j_ghost);
                face_graph.add_face(i_domain, -1);
                face_elements.add_element(e_type, face_element.nodes.data());
                j_ghost++;
            }
        }
        assert(face_elements.size() == face_graph.size());
    }

    /*Step 5*/
    vector<FaceGraphLoc> GridCreator::create_local_face_graphs(FaceGraphGlob face_graph_glob,
                                                               Utils &utils)
    {
        vector<FaceGraphLoc> face_graphs_loc;
        for (ShortIndex r{0}; r < num_part; r++)
            face_graphs_loc.emplace_back(r);

        vector<map<Index, Index>> fIDglob2loc(num_part);
        vector<map<Index, Index>> fIDloc2glob(num_part);

        vector<Index> faces_remove;
        faces_remove.reserve(face_graph_glob.size());
        /*--------------------------------------------------------------------
        Add internal faces
        --------------------------------------------------------------------*/
        for (const auto &kv : face_graph_glob.get_cellpairs())
        {
            Index fID = kv.first;
            if (face_graph_glob.is_internal_face(fID))
            {
                Index i_glob = face_graph_glob.get_i(fID);
                Index j_glob = face_graph_glob.get_j(fID);
                Index r = utils.e2r(i_glob);
                Index i_loc = utils.eIDglob2loc(i_glob);
                Index j_loc = utils.eIDglob2loc(j_glob);
                face_graphs_loc[r].add_face(i_loc, j_loc);
                Index fIDloc = face_graphs_loc[r].size() - 1;
                assert(fIDglob2loc[r].count(fID) == 0);
                assert(fIDloc2glob[r].count(fIDloc) == 0);
                fIDglob2loc[r].emplace(fID, fIDloc);
                fIDloc2glob[r].emplace(fIDloc, fID);
                faces_remove.push_back(fID);
            }
        }
        for (Index fID : faces_remove)
            face_graph_glob.remove_face(fID);

        /*--------------------------------------------------------------------
        Add partition faces
        --------------------------------------------------------------------*/

        vector<Index> next_ghost_index(num_part);
        for (ShortIndex r{0}; r < num_part; r++)
            next_ghost_index[r] = utils.n_elements_part(r);

        for (ShortIndex r_shared{0}; r_shared < num_part; r_shared++) /*Outer loop makes sure that faces of each neigbour are grouped together*/
        {
            vector<bool> r_shares_r_shared(num_part, false);
            vector<PatchInterface> patches_int_r(num_part);

            faces_remove.clear();
            for (const auto &kv : face_graph_glob.get_cellpairs())
            {
                Index fID = kv.first;
                if (face_graph_glob.is_part_face(fID))
                {
                    Index i_glob = face_graph_glob.get_i(fID);
                    Index j_glob = face_graph_glob.get_j(fID);
                    ShortIndex r_i = utils.e2r(i_glob);
                    ShortIndex r_j = utils.e2r(j_glob);

                    assert(i_glob < j_glob && r_i < r_j);

                    if (r_i == r_shared || r_j == r_shared)
                    {
                        Index i_loc = utils.eIDglob2loc(i_glob);
                        Index j_loc = utils.eIDglob2loc(j_glob);

                        face_graphs_loc[r_i].add_face(i_loc, next_ghost_index[r_i]);
                        next_ghost_index[r_i]++;
                        face_graphs_loc[r_j].add_face(j_loc, next_ghost_index[r_j]);
                        next_ghost_index[r_j]++;

                        // /*Global face should point from rank i to rank j*/
                        // face_graphs_loc[r_i].add_face(i_loc, j_loc);
                        // face_graphs_loc[r_j].add_face(j_loc, i_loc);

                        Index fIDloc_i = face_graphs_loc[r_i].size() - 1;
                        Index fIDloc_j = face_graphs_loc[r_j].size() - 1;
                        assert(fIDglob2loc[r_i].count(fID) == 0);
                        assert(fIDglob2loc[r_j].count(fID) == 0);
                        assert(fIDloc2glob[r_i].count(fIDloc_i) == 0);
                        assert(fIDloc2glob[r_j].count(fIDloc_j) == 0);
                        fIDglob2loc[r_i].emplace(fID, fIDloc_i);
                        fIDloc2glob[r_i].emplace(fIDloc_i, fID);
                        fIDglob2loc[r_j].emplace(fID, fIDloc_j);
                        fIDloc2glob[r_j].emplace(fIDloc_j, fID);

                        assert(r_shares_r_shared[r_i] == r_shares_r_shared[r_j]);
                        if (r_shares_r_shared[r_i] == false)
                        {
                            patches_int_r[r_i].FIRST_FACE = face_graphs_loc[r_i].size() - 1;
                            patches_int_r[r_j].FIRST_FACE = face_graphs_loc[r_j].size() - 1;
                            patches_int_r[r_i].rank_neighbour = r_j;
                            patches_int_r[r_j].rank_neighbour = r_i;
                            r_shares_r_shared[r_i] = true;
                            r_shares_r_shared[r_j] = true;
                            faces_remove.push_back(fID);
                        }
                    }
                }
            }
            for (Index fID : faces_remove)
                face_graph_glob.remove_face(fID);

            /*For each local graph set where the current neigbour PatchBoundary (the PatchBoundary connecting r_other ends)*/
            for (ShortIndex r{0}; r < num_part; r++)
                if (r_shares_r_shared[r])
                {
                    patches_int_r[r].N_FACES = face_graphs_loc[r].size() - patches_int_r[r].FIRST_FACE;
                    face_graphs_loc[r].add_patch_part(patches_int_r[r]);
                }
        }

        /*--------------------------------------------------------------------
        Add boundary external patches
        --------------------------------------------------------------------*/
        faces_remove.clear();
        for (const auto &p : face_graph_glob.get_patches_ext())
        {
            Index first = p.FIRST_FACE;
            Index n_faces = p.N_FACES;

            vector<bool> rank_shares_patch(num_part, false);
            assert(rank_shares_patch.back() == false); // just checking that this is actually a vector of false
            vector<PatchBoundary> patches_ext_r(num_part);

            for (Index fID{first}; fID < first + n_faces; fID++)
            {
                if (face_graph_glob.face_is_in_patch_ext(fID, p))
                {
                    Index i_glob = face_graph_glob.get_i(fID);
                    SignedIndex j_glob = face_graph_glob.get_j(fID);
                    assert(j_glob == -1); /*since all other faces are removed, j should now be a ghost cell*/
                    ShortIndex r = utils.e2r(i_glob);
                    Index i_loc = utils.eIDglob2loc(i_glob);
                    Index j_loc = -1;
                    face_graphs_loc[r].add_face(i_loc, j_loc);
                    Index fIDloc = face_graphs_loc[r].size() - 1;
                    fIDglob2loc[r].emplace(fID, fIDloc);
                    fIDloc2glob[r].emplace(fIDloc, fID);
#ifndef NDEBUG
                    faces_remove.push_back(fID);
#endif
                    /*The first face of rank r that is present in the patch is found and FIRST_FACE etc should be marked*/
                    if (!rank_shares_patch[r])
                    {
                        patches_ext_r[r].boundary_type = p.boundary_type;
                        patches_ext_r[r].FIRST_FACE = face_graphs_loc[r].size() - 1;
                        rank_shares_patch[r] = true;
                    }
                }
            }
            for (ShortIndex r{0}; r < num_part; r++)
                if (rank_shares_patch[r])
                {
                    patches_ext_r[r].N_FACES = face_graphs_loc[r].size() - patches_ext_r[r].FIRST_FACE;
                    face_graphs_loc[r].add_patch_bound(patches_ext_r[r]);
                }
        }
#ifndef NDEBUG
        for (Index fID : faces_remove)
            face_graph_glob.remove_face(fID);
        assert(face_graph_glob.size() == 0); /*Now all global faces should have been removed*/
#endif
        utils.set_fIDglob2loc(move(fIDglob2loc));
        utils.set_fIDloc2glob(move(fIDloc2glob));
        return face_graphs_loc;
    }

    /*Step 6*/
    void GridCreator::create_local_primal_grids(const PrimalGrid &primal_grid,
                                                const vector<FaceGraphLoc> &face_graphs_loc,
                                                Utils &utils,
                                                vector<unique_ptr<PrimalGrid>> &primal_grids_loc)
    {
        assert(NF_MPI::get_rank() == 0);

        for (ShortIndex r{0}; r < num_part; r++)
            primal_grids_loc.push_back(make_unique<PrimalGrid>());

        vector<map<Index, Index>> nIDglob2loc(num_part);
        vector<map<Index, Index>> nIDloc2glob(num_part);

        const vector<Vec3> &nodes = primal_grid.get_nodes();

        const Elements &e_vol = primal_grid.get_vol_elements();

        for (Index i{0}; i < e_vol.size(); i++)
        {
            Index r = utils.e2r(i);
            Elements &e_vol_r = primal_grids_loc[r]->get_vol_elements();
            vector<Vec3> &nodes_r = primal_grids_loc[r]->get_nodes();

            const Index *element = e_vol.get_element_nodes(i);
            ShortIndex num_nodes = e_vol.get_n_element_nodes(i);
            ElementType e_type = e_vol.get_element_type(i);

            /*--------------------------------------------------------------------
            Loops over all global node indices of an element. If a new node index is
            discovered, the point corresponding to the node index is added and
            a mapping from global indices are added.
            --------------------------------------------------------------------*/
            for (ShortIndex k{0}; k < num_nodes; k++)
            {
                Index nIDglob = element[k];
                if (nIDglob2loc[r].count(nIDglob) == 0)
                {
                    nodes_r.emplace_back(nodes[nIDglob]);
                    Index nIDloc = nodes_r.size() - 1;
                    nIDglob2loc[r].emplace(nIDglob, nIDloc);
                    nIDloc2glob[r].emplace(nIDloc, nIDglob);
                }
            }
            /*--------------------------------------------------------------------
            Copies an element from the global to the local domain, and renumbers
            the node indices from global to local values.
            --------------------------------------------------------------------*/
            e_vol_r.add_element_local(e_type, element, nIDglob2loc[r]);

            assert(e_vol_r.size() - 1 == utils.eIDglob2loc(i));
        }
        utils.set_nIDglob2loc(move(nIDglob2loc));
        utils.set_nIDloc2glob(move(nIDloc2glob));
        /*--------------------------------------------------------------------
        Adding the local face elements to each local primal grid
        --------------------------------------------------------------------*/
        const Elements &e_faces = primal_grid.get_face_elements();
        for (ShortIndex r{0}; r < num_part; r++)
        {
            Elements &e_faces_r = primal_grids_loc[r]->get_face_elements();
            const FaceGraphLoc &fg_loc = face_graphs_loc[r];
            for (const auto &kv : fg_loc.get_cellpairs())
            {
                Index fIDr = kv.first;
                Index fID = utils.fIDloc2glob(r, fIDr);
                ElementType e_type = e_faces.get_element_type(fID);
                const Index *e_nodes = e_faces.get_element_nodes(fID);
                e_faces_r.add_element_local(e_type, e_nodes, utils.get_nIDglob2loc(r));
            }
        }

        for (ShortIndex r{0}; r < num_part; r++)
        {
            primal_grids_loc[r]->get_nodes().shrink_to_fit();
            primal_grids_loc[r]->get_vol_elements().shrink_to_fit();
            primal_grids_loc[r]->get_face_elements().shrink_to_fit();
        }
    }

    /*Step 7*/
    void GridCreator::create_local_FV_grids(const vector<FaceGraphLoc> &face_graphs_loc,
                                            vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                            const Utils &utils,
                                            vector<unique_ptr<FV_Grid>> &FV_grids_loc)
    {
        assert(NF_MPI::get_rank() == 0);

        for (ShortIndex r{0}; r < num_part; r++)
            FV_grids_loc.push_back(make_unique<FV_Grid>());

        for (ShortIndex r{0}; r < num_part; r++)
        {
            const FaceGraphLoc &face_graph_r = face_graphs_loc[r];
            PrimalGrid &primal_grid_r = *primal_grids_loc[r];
            const Elements &e_vol_r = primal_grid_r.get_vol_elements();
            FV_Grid &FV_grid_r = *FV_grids_loc[r];
            Faces &faces_r = FV_grid_r.faces;
            Cells &cells_r = FV_grid_r.cells;

            /*Create cells*/
            Index num_ghost = face_graph_r.find_num_ghost_part() + face_graph_r.find_num_ghost_ext();
            cells_r.resize(e_vol_r.size() + num_ghost);

            /*Create faces*/
            faces_r.reserve(face_graph_r.size());

            Index next_ghost_index = e_vol_r.size() + face_graphs_loc[r].find_num_ghost_part();
            for (const auto &kv : face_graph_r.get_cellpairs())
            {
                Index i = kv.second.i;
                SignedIndex j = kv.second.j;
                assert(j >= 0 || j == -1);
                if (j == -1)
                {
                    j = next_ghost_index;
                    next_ghost_index++;
                }
                faces_r.cell_indices.push_back(geometry::Faces::Cellpair{i, static_cast<Index>(j)});
            }
            /*Set external and partition patches*/
            FV_grid_r.patches_interf = face_graph_r.get_patches_interf();
            FV_grid_r.patches_bound = face_graph_r.get_patches_ext();

            /*Reorder faces and face elements*/
            Elements &e_faces_r = primal_grid_r.get_face_elements();
            Index num_interior_faces = faces_r.size() - num_ghost;
            reorder_face_entities(num_interior_faces,
                                  // FV_grid_r.get_patches_interface(),
                                  FV_grid_r.get_patches_boundary(),
                                  faces_r,
                                  e_faces_r);
        }

        /*Finally the ghost indices of each partition needs to be reordered so that
        partition communication gets the correct access pattern.
        the faces at the partitions are already sorted by index i, j of the neigbouring
        domain sharing the patch should be sorted based on i*/
        // for (ShortIndex r{0}; r < num_part; r++)
        // {
        //     FV_Grid &FV_grid_r = *FV_grids_loc[r];
        //     Faces &faces = FV_grid_r.faces;
        //     auto &patches_interf = FV_grid_r.get_patches_interface();
        //     for (const auto &p : patches_interf)
        //     {
        //         FV_Grid &FV_grid_other = *FV_grids_loc[p.rank_neighbour];
        //         const Faces &faces_other = FV_grid_other.get_faces();
        //         const auto &patches_interf_other = FV_grid_other.get_patches_interface();
        //         /*Find the patch in other corresponding to r*/
        //         auto it = std::find_if(patches_interf_other.begin(), patches_interf_other.end(),
        //                                [r](const PatchInterface &p_other)
        //                                {
        //                                    return p_other.rank_neighbour == r;
        //                                });
        //         assert(it != patches_interf_other.end());
        //         auto &p_other = *it;
        //         assert(p.N_FACES == p_other.N_FACES);
        //         Index N_FACES = p.N_FACES;
        //         Index FIRST_FACE = p.FIRST_FACE;
        //         Index FIRST_FACE_other = p_other.FIRST_FACE;

        //         for (Index ij{FIRST_FACE_other}; ij < N_FACES; ij++)
        //         {
        //         }

        // vector<Index> j_r(N_FACES);
        // for (Index ij{0}; ij < N_FACES; ij++)
        // {
        //     j_r[ij] = faces.get_cell_j(ij + FIRST_FACE);
        //     assert(j_r[ij] >= faces.get_cell_j(FIRST_FACE) && j_r[ij] <= faces.get_cell_j(FIRST_FACE + N_FACES - 1));
        // }

        // vector<Index> i_other(N_FACES);
        // for (Index ij{0}; ij < N_FACES; ij++)
        // {
        //     i_other[ij] = faces_other.get_cell_i(ij + FIRST_FACE_other);
        //     assert(i_other[ij] >= faces_other.get_cell_i(FIRST_FACE_other) &&
        //            j_r[ij] <= faces_other.get_cell_i(FIRST_FACE_other + N_FACES - 1));
        // }

        // std::sort(j_r.begin(), j_r.end(), [&i_other](Index a, Index b)
        //           { return i_other[a] < i_other[b]; });

        // for (Index ij{0}; ij < N_FACES; ij++)
        //     faces.cell_indices[ij + FIRST_FACE].j = j_r[ij];
    }

    void GridCreator::reorder_face_entities(Index num_interior_faces,
                                            // const vector<PatchInterface> &patches_int,
                                            const vector<PatchBoundary> &patches_ext,
                                            Faces &faces,
                                            Elements &face_elements)
    {
        /*Sorts the faces based on the indices of the neighbour cells. We want them to be sorted in increasing order,
        with the owner cell i being prioritized first and the neighbour cell j second. As an example,
        the ordering {(0,1), (0,3), (0,2), (1,1)} should be changed to {(0,1), (0,2), (0,3), (1,1)}
        The way the grid is constructed this is allready achieved for the owner cell i. However j might need sorting.
        The comparison operator < for Face is defined for this purpose.
        The interior and each boundary PatchBoundary is sorted separately. The face elements are sorted accordingly*/

        Elements face_elements_to_sort;

        faces.sort_face_entities(0, num_interior_faces, face_elements, face_elements_to_sort);

        /*DONT SORT INTERFACE PATCHES! THIS SCREWS UP THE COMMUNICATION ORDERING. COULD
        ALTERNATIVELY SORT BOTH INTERFACES THE SAME WAY BASED ON min(i_ra, i_rb), BUT
        WON'T DO THIS BEFORE I SEE THAT ALL WORKS*/

        // for (const PatchInterface &p : patches_int)
        // {
        //     faces.sort_face_entities(p.FIRST_FACE, p.FIRST_FACE + p.N_FACES, face_elements, face_elements_to_sort);
        // }

        for (const PatchBoundary &p : patches_ext)
        {
            faces.sort_face_entities(p.FIRST_FACE, p.FIRST_FACE + p.N_FACES, face_elements, face_elements_to_sort);
        }
        face_elements = face_elements_to_sort;
    }

    void GridCreator::set_global_config_data(Config &config,
                                             unique_ptr<PrimalGrid> &primal_grid_glob,
                                             const FaceGraphGlob &face_graph_glob)
    {
        assert(NF_MPI::get_rank() == 0);
        /*Setting global grid metrics in config for rank 0*/
        Index n_nodes = primal_grid_glob->get_nodes().size();
        Index n_interior_cells = primal_grid_glob->get_vol_elements().size();
        Index n_exterior_faces = face_graph_glob.find_num_ghost_ext();
        Index n_interior_faces = primal_grid_glob->get_face_elements().size() - n_exterior_faces;
        Index n_connecitivty_indices = primal_grid_glob->get_vol_elements().get_e_ptr().back();
        config.set_grid_metrics_global(n_nodes, n_interior_cells, n_interior_faces, n_exterior_faces, n_connecitivty_indices);
    }

    /*Uses MPI to copy all the local PrimalGrids and FV_Grids from rank 0 to the other ranks, it
    also sets various data in config objects*/
    void GridCreator::send_recv_grids(Config &config,
                                      vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                      vector<unique_ptr<FV_Grid>> &FV_grids_loc,
                                      unique_ptr<PrimalGrid> &primal_grid,
                                      unique_ptr<FV_Grid> &FV_grid)
    {

        // if (NF_MPI::get_rank() == 0)
        // {
        //     MPI_DBG_WAIT;
        // }

        ShortIndex rank = NF_MPI::get_rank();
        ShortIndex size = NF_MPI::get_size();
        assert(size == num_part);
        if (rank == 0)
        {
            primal_grid = move(primal_grids_loc[0]);
            FV_grid = move(FV_grids_loc[0]);

            // PatchBoundary pb0;
            // pb0.boundary_type = BoundaryType::FarField;
            // pb0.FIRST_FACE = 5;
            // pb0.N_FACES = 13;
            // string bytz;
            // serialization::serialize(bytz, pb0);
            // PatchBoundary pb1;
            // serialization::deserialize(bytz, pb1);

            // string bytes;
            // serialization::serialize(bytes, *primal_grid);
            // Index sz = bytes.size();
            // PrimalGrid pg;
            // serialization::deserialize(bytes, pg);

            // string bytes;
            // serialization::serialize(bytes, *FV_grid);
            // Index sz = bytes.size();
            // FV_Grid fg;
            // serialization::deserialize(bytes, fg);

            for (ShortIndex r_loc{1}; r_loc < size; r_loc++)
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
        Index n_nodes;
        Index n_interior_cells;
        Index n_exterior_faces;
        Index n_interior_faces;
        Index n_connecitivty_indices;
        if (rank == 0)
        {
            n_nodes = config.get_N_NODES_GLOB();
            n_interior_cells = config.get_N_CELLS_INT_GLOB();
            n_exterior_faces = config.get_N_FACES_EXT_GLOB();
            n_interior_faces = config.get_N_FACES_INT_GLOB();
            n_connecitivty_indices = config.get_N_CONNECTIVITY_INDICES_GLOB();
        }
        NF_MPI::Bcast(&n_nodes, 1, 0);
        NF_MPI::Bcast(&n_interior_cells, 1, 0);
        NF_MPI::Bcast(&n_exterior_faces, 1, 0);
        NF_MPI::Bcast(&n_interior_faces, 1, 0);
        NF_MPI::Bcast(&n_connecitivty_indices, 1, 0);
        if (rank != 0)
            config.set_grid_metrics_global(n_nodes, n_interior_cells, n_interior_faces, n_exterior_faces, n_connecitivty_indices);

        NF_MPI::Barrier();

        set_config_grid_data_local(config, primal_grid, FV_grid);
    }

    void GridCreator::set_config_grid_data_local(Config &config,
                                                 unique_ptr<PrimalGrid> &primal_grid,
                                                 unique_ptr<FV_Grid> &FV_grid)
    {
        Index n_nodes = primal_grid->get_nodes().size();
        Index n_interior_cells = primal_grid->get_vol_elements().size();
        Index n_partition_faces = FV_grid->find_num_ghost_part();
        Index n_exterior_faces = FV_grid->find_num_ghost_ext();
        Index n_interior_faces = FV_grid->get_faces().size() - n_partition_faces - n_exterior_faces;
        config.set_grid_metrics_local(n_nodes, n_interior_cells, n_interior_faces, n_partition_faces, n_exterior_faces);
    }
}