#include "../../include/parallelization/GridCreator.hpp"

/*--------------------------------------------------------------------
This function reads the mesh file and creates both a primal grid
and a finite volume grid. If multiple processors are used, Metis
is used to partition the mesh.
--------------------------------------------------------------------*/

namespace geometry
{

    void GridCreator::create_global_face_entities(const Elements &vol_elements_glob,
                                                  const Vector<ElementPatch> &element_patches,
                                                  map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                                  Elements &face_elements_glob,
                                                  map<FaceElement, GhostDataPartition> &internal_boundary_faces,
                                                  const Vector<Index> &part)
    {
        if (NF_MPI::get_rank() != 0)
            return;

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
        for (const auto &element_patch : element_patches)
        {

            const Elements &surface_elements = element_patch.boundary_elements;
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
                                               Index rank_loc,
                                               const Vector<Index> &part,
                                               unique_ptr<PrimalGrid> &primal_grid_loc,
                                               map<Index, Index> &nID_glob_to_loc,
                                               map<Index, Index> &nID_loc_to_glob,
                                               map<Index, Index> &eID_glob_to_loc)
    {

        if (NF_MPI::get_rank() != 0)
            return;
        Vector<Vec3> nodes_loc;
        Elements vol_elements_loc;
        Vector<ElementPatch> element_patches_loc;

        assert(part.size() == vol_elements_glob.size());
        assert(vol_elements_loc.size() == 0);
        assert(nodes_loc.size() == 0);
        assert(element_patches_loc.size() == 0);
        assert(nID_glob_to_loc.size() == 0);

        /*Count num nodes for each rank (only used for preallocation)*/
        Index num_vol_elements_loc{0};
        for (Index i : part)
            if (i == rank_loc)
                num_vol_elements_loc++;
        vol_elements_loc.reserve(num_vol_elements_loc, MAX_NODES_VOLUME_ELEMENT);

        nodes_loc.reserve(1.5 * nodes_glob.size() / NF_MPI::get_size()); /*Assuming equal partition times a factor as estimate for number of nodes for each processor*/

        /*--------------------------------------------------------------------
        Looping over volume elements to copy from global elements and nodes
        to local elements and nodes.
        --------------------------------------------------------------------*/
        for (Index i{0}; i < vol_elements_glob.size(); i++)
        {
            if (part[i] == rank_loc)
            {
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

                assert(eID_glob_to_loc.count(i) == 0);
                eID_glob_to_loc.emplace(vol_elements_loc.size());

                vol_elements_loc.add_element_local(e_type, element, nID_glob_to_loc);
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
        vol_elements_loc.shrink_to_fit();
        nodes_loc.shrink_to_fit();
        primal_grid_loc = make_unique<PrimalGrid>(nodes_loc, vol_elements_loc, element_patches_loc);
    }

    void GridCreator::create_FV_grid_local(const Config &config,
                                           const map<FaceElement, pair<Index, long int>> &faces_to_cells_glob,
                                           const PartitionUtils &utils,
                                           Index rank_loc,
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
        Vector<Patch> patches_loc;
        Vector<PartitionPatch> part_patches_loc;

        /*--------------------------------------------------------------------
        Step 1: Create all cells from elements.
        --------------------------------------------------------------------*/

        cells_loc.resize(vol_elements_loc.size() + FV_Grid::find_N_GHOST_cells(element_patches_loc));

        /*--------------------------------------------------------------------
        Adding internal faces first
        --------------------------------------------------------------------*/
        for (const auto &p : faces_to_cells_glob)
        {
            Index u_glob = p.second.first;
            Index v_glob = p.second.second;
            const FaceElement &fe_glob = p.first;

            bool u_is_in_domain = utils.get_partID_from_global_eID(u_glob) == rank_loc;
            bool v_is_in_domain = utils.get_partID_from_global_eID(v_glob) == rank_loc;

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
        Step 4: Create the faces at the boundaries and ghost cells.
        Also add the face elements at the boundaries.
        --------------------------------------------------------------------*/
        // cout << "Creating ghost cells..\n";
        Index j_ghost = vol_elements_loc.size();
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

        /*--------------------------------------------------------------------
        Step 5: Assign interior ghost cells between domains.
        Now all faces except the ones arising from interprocessor domain
        boundaries should be accounted for. Adding each domain sequentially so
        that faces are grouped together.
        --------------------------------------------------------------------*/
        for (Index rank_neigbour{0}; rank_neigbour < NF_MPI::get_size(); rank_neigbour++)
        {
            Index num_patch_faces = 0;
            for (const auto &p : internal_boundary_faces_glob)
            {

                const GhostDataPartition &gd = p.second;

                if ((gd.rank_a == rank_loc && rank_neigbour == gd.rank_b) ||
                    (gd.rank_b == rank_loc && rank_neigbour == gd.rank_a))
                {
                    num_patch_faces++;
                    Index i_loc, j_loc;
                    assert(rank_loc == gd.rank_a || rank_loc == gd.rank_b);
                    assert(rank_neigbour == gd.rank_a || rank_neigbour == gd.rank_b);
                    if (rank_loc == gd.rank_a)
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
            PartitionPatch part_patch;
            part_patch.FIRST_FACE = faces_loc.size();
            part_patch.N_FACES = num_patch_faces;
            part_patch.rank_neighbour = rank_neigbour;
            part_patches_loc.emplace_back(part_patch);
        }
        /*-------------------------------------------------------------------
            Set some grid metrics in the config object
        --------------------------------------------------------------------*/

        faces_loc.resize_geometry_properties();

        Index num_interior_faces_loc = faces_loc.size() - FV_Grid::find_num_ghost_tot(patches_loc, part_patches_loc);
        reorder_faces_enitities(num_interior_faces_loc, patches_loc, faces_loc, face_elements_loc);

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
        assert(part_patches_loc.size() >= rank_loc && part_patches_loc.size() <= NF_MPI::get_size());

        FV_grid_loc = make_unique<FV_Grid>(cells_loc, faces_loc, patches_loc, part_patches_loc);

        // cout << "Computational grid has been created.\n";
    }

    void GridCreator::create_partitioned_grids(Config &config,
                                               unique_ptr<PrimalGrid> &primal_grid,
                                               unique_ptr<FV_Grid> &FV_grid)
    {

        const ShortIndex num_procs = NF_MPI::get_size();
        const ShortIndex rank = NF_MPI::get_rank();

        Vector<unique_ptr<PrimalGrid>> primal_grids_loc(num_procs);
        Vector<unique_ptr<FV_Grid>> FV_grids_loc(num_procs);
        if (rank == 0)
        {
            auto primal_grid_glob = make_unique<PrimalGrid>(config);

            map<FaceElement, pair<Index, long int>> faces_to_cells_glob;
            Elements face_elements_glob;

            const Vector<Index> part = NF_METIS::calc_element_partition(*primal_grid_glob, num_procs);

            map<FaceElement, GhostDataPartition> internal_ghost_faces;

            create_global_face_entities(primal_grid_glob->get_vol_elements(),
                                        primal_grid_glob->get_element_patches(),
                                        faces_to_cells_glob,
                                        face_elements_glob,
                                        internal_ghost_faces,
                                        part);

            map<Index, Index> eID_glob_to_loc;
            Vector<map<Index, Index>> nID_glob_to_loc_vec(num_procs);
            Vector<map<Index, Index>> nID_loc_to_glob_vec(num_procs);

            /*--------------------------------------------------------------------
            Build local primal grids from global primal grid (and some addressing
            structures)
            --------------------------------------------------------------------*/
            for (Index i_part{0}; i_part < num_procs; i_part++)
            {
                Elements vol_elements_loc;
                Vector<Vec3> nodes_loc;
                Vector<ElementPatch> element_patches_loc;

                create_primal_grid_local(primal_grid_glob->get_vol_elements(),
                                         primal_grid_glob->get_nodes(),
                                         primal_grid_glob->get_element_patches(),
                                         i_part,
                                         part,
                                         primal_grids_loc[i_part],
                                         nID_glob_to_loc_vec[i_part],
                                         nID_loc_to_glob_vec[i_part],
                                         eID_glob_to_loc);

                primal_grids_loc.emplace_back(nodes_loc, vol_elements_loc, element_patches_loc);
            }
            /*--------------------------------------------------------------------
            Build local FV_Grids from local primal grids
            --------------------------------------------------------------------*/

            for (Index i_part{0}; i_part < num_procs; i_part++)
            {
                PartitionUtils utils{part,
                                     eID_glob_to_loc,
                                     nID_glob_to_loc_vec[i_part], nID_loc_to_glob_vec[i_part]};
                create_FV_grid_local(config,
                                     faces_to_cells_glob,
                                     utils,
                                     i_part,
                                     FV_grids_loc[i_part],
                                     *primal_grids_loc[i_part],
                                     internal_ghost_faces);
            }
            assert(primal_grids_loc.size() == num_procs && FV_grids_loc.size() == num_procs);
        }
        scatter_grids(primal_grids_loc, FV_grids_loc, config, primal_grid, FV_grid);
    }

    void GridCreator::reorder_faces_enitities(Index num_interior_faces,
                                              const Vector<Patch> &patches,
                                              Faces &faces, Elements &face_elements)
    {
        /*Sorts the faces based on the indices of the neighbour cells. We want them to be sorted in increasing order,
        with the owner cell i being prioritized first and the neighbour cell j second. As an example,
        the ordering {(0,1), (0,3), (0,2), (1,1)} should be changed to {(0,1), (0,2), (0,3), (1,1)}
        The way the grid is constructed this is allready achieved for the owner cell i. However j might need sorting.
        The comparison operator < for Face is defined for this purpose.
        The interior and each boundary patch is sorted separately. The face elements are sorted accordingly*/

        Elements face_elements_to_sort;

        faces.sort_face_entities(0, num_interior_faces, face_elements, face_elements_to_sort);

        for (const Patch &patch : patches)
        {
            faces.sort_face_entities(patch.FIRST_FACE, patch.FIRST_FACE + patch.N_FACES, face_elements, face_elements_to_sort);
        }
        face_elements = face_elements_to_sort;
    }

    /*Uses MPI to copy all the local PrimalGrids and FV_Grids from rank 0 to the other ranks, it
    also sets various data in config objects*/
    void GridCreator::scatter_grids(Vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                    Vector<unique_ptr<FV_Grid>> &FV_grids_loc,
                                    Config &config,
                                    unique_ptr<PrimalGrid> &primal_grid,
                                    unique_ptr<FV_Grid> &FV_grid)
    {
        ShortIndex rank = NF_MPI::get_rank();
        ShortIndex num_procs = NF_MPI::get_size();

        if (rank == 0)
        {
            primal_grid = move(primal_grids_loc[0]);
            FV_grid = move(FV_grids_loc[0]);
            for (ShortIndex rank_loc{1}; rank_loc < num_procs; rank_loc++)
            {
                /*Send primal grid*/
                string bytes;
                serialization::serialize(bytes, primal_grids_loc[rank_loc]);
                Index num_bytes = bytes.size();
                NF_MPI::Send(&num_bytes, 1, rank_loc);
                NF_MPI::Send(bytes.data(), num_bytes, rank_loc);

                /*Send FV grid */
                bytes.clear();
                serialization::serialize(bytes, FV_grids_loc[rank_loc]);
                num_bytes = bytes.size();
                NF_MPI::Send(&num_bytes, 1, rank_loc);
                NF_MPI::Send(bytes.data(), num_bytes, rank_loc);
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
    }
}
