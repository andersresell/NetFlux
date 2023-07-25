#include "../../include/geometry/FV_Grid.hpp"

namespace geometry
{

    FV_Grid::FV_Grid(Config &config, PrimalGrid &primal_grid)
    {
        try
        {
            create_face_structure(config, primal_grid);
            calc_geometry_properties(config, primal_grid);
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error("Error creating grid:\n" + string(e.what()));
        }
    }
    /*Creates the finite volume face-based graph structure (each face stores indices
    to it's two neigbour cells). Additionally it saves the element nodes of each face.*/
    void FV_Grid::create_face_structure(Config &config, PrimalGrid &primal_grid)
    {
        const Vector<Vec3> &nodes = primal_grid.get_nodes();
        const Elements &volume_elements = primal_grid.get_volume_elements();
        Elements &face_elements = primal_grid.get_face_elements();
        const Vector<ElementPatch> &element_patches = primal_grid.get_element_patches();

        /*--------------------------------------------------------------------
        Step 1: Create all cells from elements.
        --------------------------------------------------------------------*/

        cells.resize(volume_elements.size() + find_N_GHOST_cells(element_patches));

        /*--------------------------------------------------------------------
        Step 2: Associate the face nodes with its two neighboring cells.
        --------------------------------------------------------------------*/

        map<SortedFaceElement, pair<Index, long int>> faces_to_cells;
        constexpr int CELL_NOT_YET_ASSIGNED{-1};

        for (Index cell_index{0}; cell_index < volume_elements.size(); cell_index++)
        {
            ElementType volume_element_type = volume_elements.get_element_type(cell_index);
            assert(is_volume_element.at(volume_element_type));
            for (ShortIndex k{0}; k < num_nodes_in_element.at(volume_element_type); k++)
            {
                FaceElement face_element = get_face_element_k_of_volume_element(volume_element_type,
                                                                                volume_elements.get_element_nodes(cell_index),
                                                                                k);
                SortedFaceElement sorted_face_element{face_element};
                assert(faces_to_cells.count(sorted_face_element) <= 1);
                if (faces_to_cells.count(sorted_face_element) == 0)
                {
                    // Discovered a new face
                    faces_to_cells.emplace(sorted_face_element, pair{cell_index, CELL_NOT_YET_ASSIGNED});
                }
                else
                {
                    assert(faces_to_cells.at(sorted_face_element).second == CELL_NOT_YET_ASSIGNED);
                    // Face allready discovered by previous cell
                    faces_to_cells.at(sorted_face_element).second = cell_index;
                }
            }
        }
        face_elements.reserve(faces_to_cells.size(), MAX_NODES_FACE_ELEMENT);
        faces.reserve(faces_to_cells.size());

        /*--------------------------------------------------------------------
        Step 3: Adding internal faces to the map. (Boundary faces are added
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
                faces.cell_indices.emplace_back(cell_i, static_cast<Index>(cell_j));
                face_elements.add_element(face.first.e_type, face.first.nodes.data());
            }
        }
        assert(face_elements.size() == faces.size());

        /*--------------------------------------------------------------------
        Step 4: Create the faces at the boundaries and ghost cells.
        Also add the face elements at the boundaries.
        --------------------------------------------------------------------*/
        cout << "Creating ghost cells..\n";
        Index j_ghost = volume_elements.size();
        for (const auto &element_patch : element_patches)
        {
            Patch p;
            p.boundary_type = config.get_boundary_type(element_patch.patch_name);
            p.FIRST_FACE = faces.size();
            p.N_FACES = element_patch.boundary_elements.size();
            patches.push_back(p);

            const Elements &surface_elements = element_patch.boundary_elements;
            for (Index ij{0}; ij < surface_elements.size(); ij++)
            {
                ElementType e_type = surface_elements.get_element_type(ij);
                FaceElement face_element{e_type, surface_elements.get_element_nodes(ij)};
                SortedFaceElement sorted_face_element{face_element};
                assert(faces_to_cells.at(face_element).second == CELL_NOT_YET_ASSIGNED);
                faces_to_cells.at(sorted_face_element).second = j_ghost;
                Index i_domain = faces_to_cells.at(sorted_face_element).first;
                assert(i_domain < j_ghost);
                faces.cell_indices.emplace_back(i_domain, j_ghost);
                face_elements.add_element(e_type, face_element.nodes.data());
                j_ghost++;
            }
        }
        assert(face_elements.size() == faces.size());

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
        print_grid(config, face_elements);

        cout << "Reorder faces..\n";
        reorder_faces(config, face_elements);

        print_grid(config, face_elements);
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

    void FV_Grid::calc_geometry_properties(const Config &config, const PrimalGrid &primal_grid)
    {
        cout << "Assigning geometrical properties to cells and faces\n";

        const Vector<Vec3> &nodes = primal_grid.get_nodes();
        const Elements &volume_elements = primal_grid.get_volume_elements();
        const Elements &face_elements = primal_grid.get_face_elements();

        assert(config.get_N_INTERIOR_CELLS() == volume_elements.size());
        assert(config.get_N_TOTAL_FACES() == face_elements.size());
        assert(config.get_N_TOTAL_CELLS() == cells.size());

        Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();

        /*Calculate properties for the cells*/
        for (Index i{0}; i < N_INTERIOR_CELLS; i++)
            volume_element_calc_geometry_properties(volume_elements.get_element_type(i),
                                                    volume_elements.get_element_nodes(i),
                                                    nodes, cells.volumes[i],
                                                    cells.centroids[i]);

        Index max_j{0}; // For consistency checking

        /*Loop over all faces to assign face normals cell-centroid-to-face-centroid vectors and ghost centroid */
        for (Index ij{0}; ij < config.get_N_TOTAL_FACES(); ij++)
        {
            Index i = faces.get_cell_i(ij);
            Index j = faces.get_cell_j(ij);
            max_j = std::max(max_j, j);
            assert(i < N_INTERIOR_CELLS); // i should never belong to a ghost cell (normal pointing from domain (i), to ghost (j))
            if (j >= N_INTERIOR_CELLS)
            {
                cells.volumes[j] = 0.0;
                calc_ghost_centroid(face_elements.get_element_type(ij),
                                    face_elements.get_element_nodes(ij),
                                    nodes,
                                    cells.centroids[i],
                                    cells.centroids[j]);
            }
            calc_face_properties(face_elements.get_element_type(ij),
                                 face_elements.get_element_nodes(ij),
                                 nodes,
                                 cells.centroids[i],
                                 cells.centroids[j],
                                 faces.face_normals[ij],
                                 faces.centroid_to_face_i[ij],
                                 faces.centroid_to_face_j[ij]);
        }
        assert(max_j == config.get_N_TOTAL_CELLS() - 1);
        /*Remember to never make an assertion like this: max_i == config.get_N_INTERIOR_CELLS() - 1.
         the max value of i may be lower than N_INTERIOR CELLS-1*/
    }

    void FV_Grid::reorder_faces(const Config &config, Elements &face_elements)
    {
        /*Sorts the faces based on the indices of the neighbour cells. We want them to be sorted in increasing order,
        with the owner cell i being prioritized first and the neighbour cell j second. As an example,
        the ordering {(0,1), (0,3), (0,2), (1,1)} should be changed to {(0,1), (0,2), (0,3), (1,1)}
        The way the grid is constructed this is allready achieved for the owner cell i. However j might need sorting.
        The comparison operator < for Face is defined for this purpose.
        The interior and each boundary patch is sorted separately. The face elements are sorted accordingly*/

        Elements face_elements_to_sort;

        Index N_INTERIOR_FACES = config.get_N_INTERIOR_FACES();

        faces.sort(0, N_INTERIOR_FACES, face_elements, face_elements_to_sort);

        for (Patch &patch : patches)
        {
            faces.sort(patch.FIRST_FACE, patch.FIRST_FACE + patch.N_FACES, face_elements, face_elements_to_sort);
        }
        face_elements = face_elements_to_sort;
    }

    Index FV_Grid::find_N_GHOST_cells(const Vector<ElementPatch> &element_patches)
    {
        Index N_GHOST{0};
        for (const auto &element_patch : element_patches)
            N_GHOST += element_patch.boundary_elements.size();
        return N_GHOST;
    }

    void FV_Grid::calc_face_properties(ElementType e_type,
                                       const Index *element,
                                       const Vector<Vec3> &nodes,
                                       const Vec3 &cell_center_i,
                                       const Vec3 &cell_center_j,
                                       Vec3 &S_ij,
                                       Vec3 &centroid_to_face_i,
                                       Vec3 &centroid_to_face_j)
    {
        face_element_calc_face_normal(e_type, element, nodes, S_ij);
        Scalar normal_dot_product = S_ij.dot(cell_center_j - cell_center_i);
        assert(normal_dot_product != 0); // Just banning this for now, altough it is possibly possible with a high skewness, but valid mesh
        if (normal_dot_product < 0)
            S_ij *= -1; // Flipping normal if it's not pointing from i to j
        Vec3 face_centroid;
        face_element_calc_centroid(e_type, element, nodes, face_centroid);
        centroid_to_face_i = face_centroid - cell_center_i;
        centroid_to_face_j = face_centroid - cell_center_j;
    }

    void FV_Grid::calc_ghost_centroid(ElementType boundary_e_type,
                                      const Index *boundary_element,
                                      const Vector<Vec3> &nodes,
                                      const Vec3 &centroid_i,
                                      Vec3 &centroid_ghost)
    {
        assert(!is_volume_element.at(boundary_e_type));
        Vec3 centroid_face;
        face_element_calc_centroid(boundary_e_type, boundary_element, nodes, centroid_face);
        centroid_ghost = 2 * centroid_face - centroid_i;
    }

    void FV_Grid::print_grid(const Config &config, const Elements &face_elements) const
    {

        cout << "CELLS:\n";
        for (Index i{0}; i < cells.size(); i++)
        {
            if (i >= config.get_N_INTERIOR_CELLS())
                cout << "GHOST, ";
            cout << cells.to_string(i) << endl;
        }

        cout << "\n\nFACES:\n";
        for (Index ij{0}; ij < faces.size(); ij++)
            cout << faces.to_string(ij) << ", Element: " << face_elements.to_string(ij) << endl;

        cout << "\n\nPATCHES:\n";
        for (const auto &patch : patches)
        {
            cout << "patch type: " << (int)patch.boundary_type << "\nFIRST FACE: " << patch.FIRST_FACE << "\nN_FACES: " << patch.N_FACES << "\n\n";
        }
    }

}