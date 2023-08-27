#include "../../include/geometry/FV_Grid.hpp"

namespace geometry
{

    // FV_Grid::FV_Grid(Cells &&cells, Faces &&faces, vector<PatchBoundary> &&patches_bound, vector<PatchInterface> &&partition_patchespatches_interf)

    //     : cells{cells}, faces{faces}, patches_bound{patches_bound}, patches_interf{patches_interf}
    // {
    // }

    // FV_Grid::FV_Grid(Config &config, PrimalGrid &primal_grid)
    // {
    //     try
    //     {
    //         create_face_structure(config, primal_grid);
    //         calc_geometry_properties(config, primal_grid);
    //         // print_grid(config, primal_grid.face_elements);
    //     }
    //     catch (const std::exception &e)
    //     {
    //         throw std::runtime_error("Error creating grid:\n" + string(e.what()));
    //     }
    // }

    /*Creates the finite volume face-based graph structure (each face stores indices
    to it's two neigbour cells). Additionally it saves the element nodes of each face.*/
    // void FV_Grid::create_face_structure(Config &config, PrimalGrid &primal_grid)
    // {
    //     const vector<Vec3> &nodes = primal_grid.get_nodes();
    //     const Elements &vol_elements = primal_grid.get_vol_elements();
    //     Elements &face_elements = primal_grid.get_face_elements();
    //     const vector<ElementPatch> &element_patches = primal_grid.get_element_patches();

    //     /*--------------------------------------------------------------------
    //     Step 1: Create all cells from elements.
    //     --------------------------------------------------------------------*/

    //     cells.resize(vol_elements.size() + find_N_GHOST_cells(element_patches));

    //     /*--------------------------------------------------------------------
    //     Step 2: Associate the face nodes with its two neighboring cells.
    //     --------------------------------------------------------------------*/

    //     map<FaceElement, pair<Index, long int>> faces_to_cells;
    //     constexpr int CELL_NOT_YET_ASSIGNED{-1};

    //     for (Index cell_index{0}; cell_index < vol_elements.size(); cell_index++)
    //     {
    //         ElementType volume_element_type = vol_elements.get_element_type(cell_index);
    //         assert(is_volume_element(volume_element_type));
    //         for (ShortIndex k{0}; k < get_num_faces_volume_element(volume_element_type); k++)
    //         {
    //             FaceElement face_element = get_face_element_k_of_volume_element(volume_element_type,
    //                                                                             vol_elements.get_element_nodes(cell_index),
    //                                                                             k);

    //             assert(faces_to_cells.count(face_element) <= 1);
    //             if (faces_to_cells.count(face_element) == 0)
    //             {
    //                 // Discovered a new face
    //                 faces_to_cells.emplace(face_element, pair{cell_index, CELL_NOT_YET_ASSIGNED});
    //             }
    //             else
    //             {
    //                 assert(faces_to_cells.at(face_element).second == CELL_NOT_YET_ASSIGNED);
    //                 // Face allready discovered by previous cell
    //                 faces_to_cells.at(face_element).second = cell_index;
    //             }
    //         }
    //     }
    //     face_elements.reserve(faces_to_cells.size(), MAX_NODES_FACE_ELEMENT);
    //     faces.reserve(faces_to_cells.size());

    //     /*--------------------------------------------------------------------
    //     Step 3: Adding internal faces to the map. (Boundary faces are added
    //     later, this is to get the correct grouping of patches). Face elements
    //     (used for calculating geometry properties) are created simultaneously
    //     as faces, to get the correct ordering.
    //     --------------------------------------------------------------------*/

    //     for (const auto &face : faces_to_cells)
    //     {
    //         Index cell_i = face.second.first;
    //         long int cell_j = face.second.second;
    //         if (cell_j != CELL_NOT_YET_ASSIGNED)
    //         { // Only add internal faces
    //             assert(cell_i < cell_j);
    //             faces.cell_indices.emplace_back(cell_i, static_cast<Index>(cell_j));
    //             face_elements.add_element(face.first.e_type, face.first.nodes.data());
    //         }
    //     }
    //     assert(face_elements.size() == faces.size());

    //     /*--------------------------------------------------------------------
    //     Step 4: Create the faces at the boundaries and ghost cells.
    //     Also add the face elements at the boundaries.
    //     --------------------------------------------------------------------*/
    //     cout << "Creating ghost cells..\n";
    //     Index j_ghost = vol_elements.size();
    //     for (const auto &element_PatchBoundary : element_patches)
    //     {
    //         PatchBoundary p;
    //         p.boundary_type = config.get_boundary_type(element_PatchBoundary.PatchBoundary_name);
    //         p.FIRST_FACE = faces.size();
    //         p.N_FACES = element_PatchBoundary.boundary_elements.size();
    //         patches.push_back(p);

    //         const Elements &surface_elements = element_PatchBoundary.boundary_elements;
    //         for (Index ij{0}; ij < surface_elements.size(); ij++)
    //         {
    //             ElementType e_type = surface_elements.get_element_type(ij);
    //             FaceElement face_element{e_type, surface_elements.get_element_nodes(ij)};
    //             assert(faces_to_cells.at(face_element).second == CELL_NOT_YET_ASSIGNED);
    //             faces_to_cells.at(face_element).second = j_ghost;
    //             Index i_domain = faces_to_cells.at(face_element).first;
    //             assert(i_domain < j_ghost);
    //             faces.cell_indices.emplace_back(i_domain, j_ghost);
    //             face_elements.add_element(e_type, face_element.nodes.data());
    //             j_ghost++;
    //         }
    //     }
    //     assert(face_elements.size() == faces.size());

    //     /*--------------------------------------------------------------------
    //     Step 5: Assign interior ghost cells between domains.
    //     Now all faces except the ones arising from interprocessor domain
    //     boundaries should be accounted for.
    //     --------------------------------------------------------------------*/

    //     /*-------------------------------------------------------------------
    //         Set some grid metrics in the config object
    //     --------------------------------------------------------------------*/
    //     Index N_NODES = nodes.size();
    //     Index N_INTERIOR_CELLS = vol_elements.size();
    //     Index N_TOTAL_CELLS = cells.size();
    //     Index N_INTERIOR_FACES = faces.size() - find_N_GHOST_cells(element_patches);
    //     Index N_TOTAL_FACES = faces.size();
    //     config.set_grid_metrics(N_NODES, N_INTERIOR_CELLS, N_TOTAL_CELLS, N_INTERIOR_FACES, N_TOTAL_FACES);

    //     /*--------------------------------------------------------------------
    //     Sort faces so that the interior faces appear first with the owner index
    //     i allways being less than neigbour index j. The same logic is applied
    //     PatchBoundary-wise to the boundaries
    //     --------------------------------------------------------------------*/

    //     /*All Vectors within faces and cells must have the correct size before reordering*/
    //     faces.resize_geometry_properties();

    //     reorder_faces(config, face_elements);

    //     /*--------------------------------------------------------------------
    //     Ensuring that no unneccessary memory isn't used
    //     --------------------------------------------------------------------*/
    //     primal_grid.element_patches.clear(); // element patches no longer needed.
    //     assert(primal_grid.element_patches.empty());
    //     assert(faces.cell_indices.capacity() == N_TOTAL_FACES);
    //     assert(faces.face_normals.capacity() == N_TOTAL_FACES);
    //     assert(faces.centroid_to_face_i.capacity() == N_TOTAL_FACES);
    //     assert(faces.centroid_to_face_j.capacity() == N_TOTAL_FACES);
    //     assert(cells.volumes.capacity() == N_TOTAL_CELLS);
    //     assert(cells.centroids.capacity() == N_TOTAL_CELLS);

    //     cout << "Computational grid has been created.\n";
    // }

    void FV_Grid::calc_geometry_properties(const Config &config, const PrimalGrid &primal_grid, PartitionComm &part_comm)
    {
        cout << "Assigning geometrical properties to cells and faces.\n";

        const vector<Vec3> &nodes = primal_grid.get_nodes();
        const Elements &vol_elements = primal_grid.get_vol_elements();
        const Elements &face_elements = primal_grid.get_face_elements();

        assert(config.get_N_CELLS_INT() == vol_elements.size());
        assert(config.get_N_FACES_TOT() == face_elements.size());
        assert(config.get_N_CELLS_TOT() == cells.size());

        /*NEED TO UPDATE THIS TO ALSO CALCULATE CENTROID OF PARTITION GHOST CELL.
        WILL POSSIBLY NEED TO EDIT THE CODE UPSTREAM SO THAT THE GHOST ELEMENT IS COPIED,
        SINCE THIS INFO IS NEEDED TO CALCULATE THE CENTROID.*/
        Index N_CELLS_INT = config.get_N_CELLS_INT();
        Index first_face_interf;
        bool patches_interf_exist = false;
        if (!patches_interf.empty())
        {
            first_face_interf = patches_interf[0].FIRST_FACE;
            patches_interf_exist = true;
        }
        Index first_face_bound;
        bool patches_bound_exist = false;
        if (!patches_bound.empty())
        {
            first_face_bound = patches_bound[0].FIRST_FACE;
            patches_bound_exist = true;
        }

#ifndef NDEBUG
        if (patches_interf_exist && patches_bound_exist)
            assert(first_face_bound > first_face_interf);
#endif
        /*Calculate properties for the interior cells*/
        for (Index i{0}; i < N_CELLS_INT; i++)
            volume_element_calc_geometry_properties(vol_elements.get_element_type(i),
                                                    vol_elements.get_element_nodes(i),
                                                    nodes, cells.volumes[i],
                                                    cells.centroids[i]);

        Index max_j{0}; // For consistency checking

        /*Loop over all faces to assign face normals cell-centroid-to-face-centroid vectors and ghost centroid */
        for (Index ij{0}; ij < config.get_N_FACES_TOT(); ij++)
        {
            Index i = faces.get_cell_i(ij);
            Index j = faces.get_cell_j(ij);
            assert(i != j);
            max_j = std::max(max_j, j);
            assert(i < N_CELLS_INT); // i should never belong to a ghost cell (normal pointing from domain (i), to ghost (j))
            if (patches_interf_exist && (ij >= first_face_interf) && (ij < first_face_bound))
            {
                cells.volumes[j] = 0.0;
            }
            else if (patches_bound_exist && ij >= first_face_bound)
            {
                assert(j >= config.get_N_CELLS_DOMAIN() && j < config.get_N_CELLS_TOT());
                cells.volumes[j] = 0.0;
                calc_ghost_centroid_bound(face_elements.get_element_type(ij),
                                          face_elements.get_element_nodes(ij),
                                          nodes,
                                          cells.centroids[i],
                                          cells.centroids[j]);
            }
        }
        assert(max_j == config.get_N_CELLS_TOT() - 1);

        part_comm.communicate_interface_ghost_centroids(cells.centroids);

        for (Index ij{0}; ij < config.get_N_FACES_TOT(); ij++)
        {
            Index i = faces.get_cell_i(ij);
            Index j = faces.get_cell_j(ij);
            calc_face_properties(face_elements.get_element_type(ij),
                                 face_elements.get_element_nodes(ij),
                                 nodes,
                                 cells.centroids[i],
                                 cells.centroids[j],
                                 faces.face_normals[ij],
                                 faces.centroid_to_face_i[ij],
                                 faces.centroid_to_face_j[ij]);
        }
    }

    Index FV_Grid::find_num_ghost_ext() const
    {
        Index N_GHOST{0};
        for (const auto &p : patches_bound)
            N_GHOST += p.N_FACES;
        assert(patches_bound.size() == 0 || N_GHOST > 0);
        return N_GHOST;
    }
    Index FV_Grid::find_num_ghost_part() const
    {
        Index N_GHOST{0};
        for (const auto &p : patches_interf)
            N_GHOST += p.N_FACES;
        assert(patches_interf.size() == 0 || N_GHOST > 0);
        return N_GHOST;
    }

    void FV_Grid::calc_face_properties(ElementType e_type,
                                       const Index *element,
                                       const vector<Vec3> &nodes,
                                       const Vec3 &cell_center_i,
                                       const Vec3 &cell_center_j,
                                       Vec3 &S_ij,
                                       Vec3 &centroid_to_face_i,
                                       Vec3 &centroid_to_face_j)
    {
        face_element_calc_face_normal(e_type, element, nodes, S_ij);
        Scalar normal_dot_product = S_ij.dot(cell_center_j - cell_center_i);
        assert(normal_dot_product != 0); // Just banning this for now, altough it might be possibly possible with a high skewness, but valid mesh
        if (normal_dot_product < 0)
            S_ij *= -1; // Flipping normal if it's not pointing from i to j
        Vec3 face_centroid;
        face_element_calc_centroid(e_type, element, nodes, face_centroid);
        centroid_to_face_i = face_centroid - cell_center_i;
        centroid_to_face_j = face_centroid - cell_center_j;
    }

    void FV_Grid::calc_ghost_centroid_bound(ElementType boundary_e_type,
                                            const Index *boundary_element,
                                            const vector<Vec3> &nodes,
                                            const Vec3 &centroid_i,
                                            Vec3 &centroid_ghost)
    {
        assert(!is_volume_element(boundary_e_type));
        Vec3 centroid_face;
        face_element_calc_centroid(boundary_e_type, boundary_element, nodes, centroid_face);
        centroid_ghost = 2 * centroid_face - centroid_i;
    }

    void FV_Grid::print_grid(const Config &config, const Elements &face_elements) const
    {

        cout << "CELLS:\n";
        for (Index i{0}; i < cells.size(); i++)
        {
            if (i >= config.get_N_CELLS_INT())
                cout << "GHOST, ";
            cout << cells.to_string(i) << endl;
        }

        cout << "\n\nFACES:\n";
        for (Index ij{0}; ij < faces.size(); ij++)
            cout << faces.to_string(ij) << ", Element: " << face_elements.to_string(ij) << endl;

        cout << "\n\npatches:\n";
        for (const auto &p : patches_bound)
        {
            cout << "PatchBoundary type: " << (int)p.boundary_type << "\nFIRST FACE: " << p.FIRST_FACE << "\nN_FACES: " << p.N_FACES << "\n\n";
        }
    }
}