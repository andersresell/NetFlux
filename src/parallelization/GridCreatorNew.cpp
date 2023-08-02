#include "../../include/parallelization/GridCreatorNew.hpp"

namespace geometry
{ /*Main loop*/
    void GridCreatorNew::create_partitioned_grids(Config &config,
                                                  unique_ptr<PrimalGrid> &primal_grid_glob,
                                                  unique_ptr<PrimalGrid> &primal_grid,
                                                  unique_ptr<FV_Grid> &FV_grid)
    {
        const ShortIndex num_procs = NF_MPI::get_size();
        const ShortIndex rank = NF_MPI::get_rank();

        if (rank == 0)
        {
            /*Step 1*/
            primal_grid_glob = make_unique<PrimalGrid>(config);
            /*Step 2*/
            const Vector<ShortIndex> part = NF_METIS::calc_element_partition(*primal_grid_glob, num_procs);
            Utils utils{part};

            /*Step 3*/
            reorder_global_grid(*primal_grid_glob, utils);

            /*Step 4*/
            FaceGraphGlob face_graph_glob{part};
            Elements face_elements_glob;
            create_global_face_graph_and_face_elements(config,
                                                       *primal_grid_glob,
                                                       face_graph_glob,
                                                       face_elements_glob,
                                                       utils);
        }
    }
    /*Step 3*/
    void GridCreatorNew::reorder_global_grid(PrimalGrid &primal_grid,
                                             Utils &utils)
    {
        assert(NF_MPI::get_rank() == 0);
        ShortIndex num_procs = NF_MPI::get_size();
        Vector<pair<Index, Index>> part2e_range(num_procs);

        Elements &vol_elements_old = primal_grid.get_vol_elements();
        Index n_elem = vol_elements_old.size();
        Vector<Index> eIDglob2loc(n_elem);
        Elements vol_elements_new;
        vol_elements_new.reserve(n_elem, MAX_NODES_VOLUME_ELEMENT);
        for (ShortIndex r_loc{0}; r_loc < num_procs; r_loc++)
        {
            /*This loop can be made faster (outermost loop removed) by first creating part_to_element_range
        (count occurences of each rank) and then using the intervals to insert elements in vol_elements_new*/

            bool first_element_found = false;
            for (Index i{0}; i < n_elem; i++)
            {
                if (utils.e2r(i) == r_loc)
                {
                    if (!first_element_found)
                        part2e_range[r_loc].first = vol_elements_new.size();
                    first_element_found = true;
                    vol_elements_new.add_element(vol_elements_old.get_element_type(i), vol_elements_old.get_element_nodes(i));

                    eIDglob2loc[i] = (r_loc == 0) ? vol_elements_new.size() : vol_elements_new.size() - part2e_range[r_loc].first;
                }
                part2e_range[r_loc].second = vol_elements_new.size();
            }
        }
        assert(vol_elements_new.size() == vol_elements_old.size());
        assert(eIDglob2loc.size() == vol_elements_old.size());
        utils.set_eIDglob2loc(move(eIDglob2loc));
        vol_elements_new.shrink_to_fit();
        vol_elements_old = vol_elements_new;
    }

    /*Step 4*/
    void GridCreatorNew::create_global_face_graph_and_face_elements(const Config &config,
                                                                    PrimalGrid &primal_grid,
                                                                    FaceGraphGlob &face_graph,
                                                                    Utils &utils)
    {
        Elements &face_elements = primal_grid.get_face_elements();
        assert(face_graph.size() == 0);
        assert(face_elements.size() == 0);

        const Vector<Vec3> &nodes = primal_grid.get_nodes();
        const Elements &vol_elements = primal_grid.get_vol_elements();
        const Vector<ElementPatch> &element_patches = primal_grid.get_element_patches();

        /*--------------------------------------------------------------------
        Associate the face nodes with its two neighboring cells.
        --------------------------------------------------------------------*/

        map<FaceElement, pair<Index, long int>> faces_to_cells;
        constexpr int CELL_NOT_YET_ASSIGNED{-1};

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
        for (const auto &element_patch : element_patches)
        {
            Patch p;
            p.boundary_type = config.get_boundary_type(element_patch.patch_name);
            p.FIRST_FACE = face_graph.size();
            p.N_FACES = element_patch.boundary_elements.size();
            face_graph.add_patch(p);

            const Elements &surface_elements = element_patch.boundary_elements;
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
    Vector<FaceGraphLoc> GridCreatorNew::create_local_face_graphs(FaceGraphGlob face_graph_glob,
                                                                  Utils &utils)
    {
        const ShortIndex size = NF_MPI::get_size();
        Vector<FaceGraphLoc> face_graphs_loc;
        for (ShortIndex r{0}; r < size; r++)
            face_graphs_loc.emplace_back(r);

        Vector<map<Index, Index>> fIDglob2loc(size);
        Vector<map<Index, Index>> fIDloc2glob(size);

        /*--------------------------------------------------------------------
        Add internal faces
        --------------------------------------------------------------------*/
        for (const auto &kv : face_graph_glob.get_cellpairs())
        {
            Index fID = kv.first;
            if (face_graph_glob.is_internal_face(fID))
            {
                Index i = face_graph_glob.get_i(fID);
                Index j = face_graph_glob.get_j(fID);
                Index r = utils.e2r(i);
                face_graphs_loc[r].add_face(i, j);
                Index loc_size = face_graphs_loc[r].size();
                assert(fIDglob2loc[r].count(fID) == 0);
                assert(fIDloc2glob[r].count(loc_size - 1) == 0);
                fIDglob2loc[r].emplace(fID, loc_size - 1);
                fIDloc2glob[r].emplace(loc_size - 1, fID);
                face_graph_glob.remove_face(fID);
            }
        }

        /*--------------------------------------------------------------------
        Add partition faces
        --------------------------------------------------------------------*/
        for (ShortIndex r_other{0}; r_other < size; r_other++) /*Outer loop makes sure that faces of each neigbour are grouped together*/
        {
            /*For each local graph set where the current neigbour patch (the patch connecting r_other begins)*/
            for (ShortIndex r{0}; r < size; r++)
                face_graphs_loc[r].set_part_begin(r_other, face_graphs_loc[r].size());

            for (const auto &kv : face_graph_glob.get_cellpairs())
            {
                Index fID = kv.first;
                if (face_graph_glob.is_part_face(fID))
                {
                    Index i = face_graph_glob.get_i(fID);
                    Index j = face_graph_glob.get_j(fID);
                    ShortIndex r_i = utils.e2r(i);
                    ShortIndex r_j = utils.e2r(j);

                    assert(i < j && r_i < r_j);

                    if (r_i == r_other || r_j == r_other)
                    {
                        /*Global face should point from rank i to rank j*/
                        face_graphs_loc[r_i].add_face(i, j);
                        face_graphs_loc[r_j].add_face(j, i);

                        Index loc_size_i = face_graphs_loc[r_i].size();
                        Index loc_size_j = face_graphs_loc[r_j].size();
                        assert(fIDglob2loc[r_i].count(fID) == 0);
                        assert(fIDglob2loc[r_j].count(fID) == 0);
                        assert(fIDloc2glob[r_i].count(loc_size_i - 1) == 0);
                        assert(fIDloc2glob[r_j].count(loc_size_j - 1) == 0);
                        fIDglob2loc[r_i].emplace(fID, loc_size_i - 1);
                        fIDloc2glob[r_i].emplace(loc_size_i - 1, fID);
                        fIDglob2loc[r_j].emplace(fID, loc_size_j - 1);
                        fIDloc2glob[r_j].emplace(loc_size_j - 1, fID);

                        face_graph_glob.remove_face(fID);
                    }
                }
            }
            /*For each local graph set where the current neigbour patch (the patch connecting r_other ends)*/
            for (ShortIndex r{0}; r < size; r++)
                face_graphs_loc[r].set_part_end(r_other, face_graphs_loc[r].size());
        }

        /*--------------------------------------------------------------------
        Add boundary patches
        --------------------------------------------------------------------*/
        for (const auto &p : face_graph_glob.get_patches())
        {
            Index first = p.FIRST_FACE;
            Index n_faces = p.N_FACES;
            BoundaryType bt = p.boundary_type;

            Vector<bool> p_already_in_r(size, false);
            assert(p_already_in_r.back() == false); // just checking that this is actually a vector of false

            for (Index fID{first}; fID < n_faces; fID++)
            {
                if (face_graph_glob.face_is_in_patch(fID, p))
                {
                    Index i = face_graph_glob.get_i(fID);
                    Index j = face_graph_glob.get_j(fID);
                    assert(j == -1); /*since all other faces are removed, j should now be a ghost cell*/
                    ShortIndex r = utils.e2r(i);
                    face_graphs_loc[r].add_face(i, j);
                    face_graph_glob.remove_face(fID);
                    if (!p_already_in_r[r])
                    {
                        face_graphs_loc[r].add_patch(bt, face_graphs_loc[r].size() - 1);
                        p_already_in_r[r] = true;
                    }
                }
            }
            for (ShortIndex r{0}; r < size; r++)
                if (p_already_in_r[r])
                    face_graphs_loc[r].set_patch_end(face_graphs_loc[r].size());
        }
        assert(face_graph_glob.size() == 0); /*Now all global faces should have been removed*/
    }

    /*Step 6*/
    Vector<PrimalGrid> GridCreatorNew::create_local_primal_grids(const PrimalGrid &primal_grid,
                                                                 const Vector<FaceGraphLoc> &face_graphs_loc,
                                                                 Utils &utils)
    {
        assert(NF_MPI::get_rank() == 0);
        const ShortIndex size = NF_MPI::get_size();
        Vector<PrimalGrid> primal_grids_loc(size);

        Vector<map<Index, Index>> nIDglob2loc(size);
        Vector<map<Index, Index>> nIDloc2glob(size);

        const Vector<Vec3> &nodes = primal_grid.get_nodes();

        const Elements &e_vol = primal_grid.get_vol_elements();

        for (Index i{0}; i < e_vol.size(); i++)
        {
            Index r = utils.e2r(i);
            Elements &e_vol_r = primal_grids_loc[r].get_vol_elements();
            Vector<Vec3> &nodes_r = primal_grids_loc[r].get_nodes();

            const Index *element = e_vol.get_element_nodes(i);
            ShortIndex num_nodes = e_vol.get_n_element_nodes(i);
            ElementType e_type = e_vol.get_element_type(i);

            /*--------------------------------------------------------------------
            Loops over all global node indices of an element. If a new node index is
            discovered, the point corresponding to the node index is added and the
            a mapping from global indices are added.
            --------------------------------------------------------------------*/
            for (ShortIndex k{0}; k < num_nodes; k++)
            {
                Index nIDglob = element[k];
                if (nIDglob2loc[r].count(nIDglob) == 0)
                {
                    nodes_r.emplace_back(nodes[nIDglob]);
                    Index nIDloc = nodes_r.size();
                    nIDglob2loc[r].emplace(nIDglob, nIDloc);
                    nIDloc2glob[r].emplace(nIDloc, nIDglob);
                }
            }
            /*--------------------------------------------------------------------
            Copies an element from the global to the local domain, and renumbers
            the node indices from global to local values.
            --------------------------------------------------------------------*/
            e_vol_r.add_element_local(e_type, element, nIDglob2loc[r]);

            assert(e_vol_r.size() == utils.eIDglob2loc(i));
        }

        /*--------------------------------------------------------------------
        Adding the local face elements to each local primal grid
        --------------------------------------------------------------------*/
        const Elements &e_faces = primal_grid.get_face_elements();
        for (ShortIndex r{0}; r < size; r++)
        {
            Elements &e_faces_r = primal_grids_loc[r].get_face_elements();
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

        for (ShortIndex r{0}; r < size; r++)
        {
            primal_grids_loc[r].get_nodes().shrink_to_fit();
            primal_grids_loc[r].get_vol_elements().shrink_to_fit();
            primal_grids_loc[r].get_face_elements().shrink_to_fit();
        }
    }