#include "../../include/geometry/FV_Grid.hpp"

using namespace geometry;

FV_Grid::FV_Grid(Config &config, const PrimalGrid &primal_grid)
{

    try
    {
        create_grid(config, primal_grid);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Error creating grid:\n" + string(e.what()));
    }
}

void FV_Grid::create_grid(Config &config, const PrimalGrid &primal_grid)
{
    cout << "Creating grid from native mesh\n";

    /*--------------------------------------------------------------------
    Step 1: Create all cells from elements, this is straight forward
    --------------------------------------------------------------------*/

    cells.reserve(tet_connect.size() + find_N_GHOST_cells());

        for (Index i{0}; i < tet_connect.size(); i++)
    {
        TetConnect tc_i = tet_connect[i];
        Tetrahedron tet = tet_from_connect(tc_i);
        cells.cell_volumes.emplace_back(tet.calc_volume());
        cells.centroids.emplace_back(tet.calc_centroid());
    }

    /*--------------------------------------------------------------------
    Step 2: Associate the face nodes with its two neighboring cells.
    --------------------------------------------------------------------*/

    Vector<Index> face_indices;
    map<SortedTriConnect, pair<Index, long int>> face_to_cells;
    constexpr int CELL_NOT_YET_ASSIGNED{-1};

    for (Index cell_index{0}; cell_index < tet_connect.size(); cell_index++)
    {

        for (ShortIndex k{0}; k < N_TET_FACES; k++)
        {
            TriConnect tc{tet_face_connectivity(tet_connect[cell_index], k)};

            // Get the sorted node indices of face k
            SortedTriConnect face_ij{tc};
            assert(face_to_cells.count(face_ij) <= 1);
            if (face_to_cells.count(face_ij) == 0)
            {
                // Discovered a new face
                face_to_cells.emplace(face_ij, pair{cell_index, CELL_NOT_YET_ASSIGNED});
            }
            else
            {

                assert(face_to_cells.at(face_ij).second == CELL_NOT_YET_ASSIGNED);

                // Face allready discovered by previous cell
                face_to_cells.at(face_ij).second = cell_index;
            }
        }
    }

    Vector<Triangle> face_triangles;
    face_triangles.reserve(face_to_cells.size());
    faces.reserve(face_to_cells.size());

    /*--------------------------------------------------------------------
    Step 3: Adding internal faces to the map. (Boundary faces are added
    later, this is to get the correct grouping of patches). Face triangles
    (used for calculating geometry properties) are created simultaneously
    as faces, to get the correct ordering.
    --------------------------------------------------------------------*/

    for (const auto &face : face_to_cells)
    {
        Index cell_i = face.second.first;
        long int cell_j = face.second.second;
        if (cell_j != CELL_NOT_YET_ASSIGNED)
        { // Only add internal faces
            assert(cell_i < cell_j);
            faces.cell_indices.emplace_back(cell_i, static_cast<Index>(cell_j));
            face_triangles.emplace_back(tri_from_connect(face.first));
        }
    }
    assert(face_triangles.size() == faces.size());

    /*--------------------------------------------------------------------
    Step 4: Create the faces at the boundaries and ghost cells.
    Also add the face triangles at the boundaries.
    --------------------------------------------------------------------*/
    cout << "Creating ghost cells..\n";
    for (const auto &tpc : tri_patch_connect_list)
    {

        Patch p;
        p.boundary_type = config.get_boundary_type(tpc.patch_name);
        p.FIRST_FACE = faces.size();
        p.N_FACES = tpc.triangles.size();
        patches.push_back(p);

        for (const TriConnect &tc : tpc.triangles)
        {
            SortedTriConnect face_ij{tc};
            Index cell_j_ghost = cells.size();
            assert(face_to_cells.at(face_ij).second == CELL_NOT_YET_ASSIGNED);
            face_to_cells.at(face_ij).second = cell_j_ghost; // this won't be used, so strictly unnecessary
            Index cell_i_domain = face_to_cells.at(face_ij).first;
            assert(cell_i_domain < cell_j_ghost);
            assert(cell_j_ghost >= tet_connect.size());
            faces.cell_indices.emplace_back(cell_i_domain, cell_j_ghost);

            face_triangles.emplace_back(tri_from_connect(tc));
            cells.add_empty(); // Adding ghost cell
        }
    }
    assert(face_triangles.size() == faces.size());

    /*-------------------------------------------------------------------
        Set some grid metrics in the config object
    --------------------------------------------------------------------*/
    Index N_NODES = nodes.size();
    Index N_INTERIOR_CELLS = tet_connect.size();
    Index N_TOTAL_CELLS = cells.size();
    Index N_INTERIOR_FACES = faces.size() - find_N_GHOST_cells();
    Index N_TOTAL_FACES = faces.size();
    config.set_grid_metrics(N_NODES, N_INTERIOR_CELLS, N_TOTAL_CELLS, N_INTERIOR_FACES, N_TOTAL_FACES);

    /*--------------------------------------------------------------------
    Assigning geometrical properties to faces and ghost cells
    --------------------------------------------------------------------*/
    cout << "Assign geometrical properties..\n";
    assign_geometry_properties(config, face_triangles);

    /*--------------------------------------------------------------------
    Sort faces so that the interior faces appear first with the owner index
    i allways being less than neigbour index j. The same logic is applied
    patch-wise to the boundaries
    --------------------------------------------------------------------*/
    cout << "Reorder faces..\n";
    reorder_faces(config);

    /*--------------------------------------------------------------------
    Reducing allocated memory
    --------------------------------------------------------------------*/
    shrink_vectors();

    cout << "Computational grid has been created.\n";
}

void FV_Grid::assign_geometry_properties(const Config &config, const Vector<Triangle> &face_triangles)
{
    faces.resize_geometry_properties();

    Index max_j{0}; // For consistency checking

    // loop over all faces to assign area vector and centroid vectors
    assert(face_triangles.size() == config.get_N_TOTAL_FACES());

    Index N_INTERIOR_CELLS = config.get_N_INTERIOR_CELLS();

    for (Index ij{0}; ij < config.get_N_TOTAL_FACES(); ij++)
    {
        Index i = faces.get_cell_i(ij);
        Index j = faces.get_cell_j(ij);
        max_j = std::max(max_j, j);

        assert(i < N_INTERIOR_CELLS); // i should never belong to a ghost cell (normal pointing outwards)

        if (j >= N_INTERIOR_CELLS)
        {
            cells.cell_volumes[j] = 0.0;
            cells.centroids[j] = calc_ghost_centroid(cells.centroids[i], face_triangles[ij]);
        }

        assign_face_properties(faces.normal_areas[ij],
                               faces.centroid_to_face_i[ij],
                               faces.centroid_to_face_j[ij],
                               face_triangles[ij],
                               cells.centroids[i],
                               cells.centroids[j]);
    }
    /*Remember to not make an assertion like this: max_i == config.get_N_INTERIOR_CELLS() - 1.
     the max value of i may be lower than N_INTERIOR CELLS-1*/

    assert(max_j == config.get_N_TOTAL_CELLS() - 1);
    assert(config.get_N_TOTAL_CELLS() == cells.size());
}

void FV_Grid::reorder_faces(const Config &config)
{
    /*Sorts the faces based on the indices of the neighbour cells. We want them to be sorted in increasing order,
    with the owner cell i being prioritized first and the neighbour cell j second. As an example,
    the ordering {(0,1), (0,3), (0,2), (1,1)} should be changed to {(0,1), (0,2), (0,3), (1,1)}
    The way the grid is constructed this is allready achieved for the owner cell i. However j might need sorting.
    The comparison operator < for Face is defined for this purpose.
    The interior and each boundary patch is sorted separately.*/
    Index N_INTERIOR_FACES = config.get_N_INTERIOR_FACES();

    faces.sort(0, N_INTERIOR_FACES);

    for (Patch &patch : patches)
    {
        faces.sort(patch.FIRST_FACE, patch.FIRST_FACE + patch.N_FACES);
    }
}

Tetrahedron FV_Grid::tet_from_connect(const TetConnect &tc) const
{
    return Tetrahedron(nodes.at(tc.a()), nodes.at(tc.b()), nodes.at(tc.c()), nodes.at(tc.d()));
}

Triangle FV_Grid::tri_from_connect(const TriConnect &tc) const
{
    return Triangle(nodes.at(tc.a()), nodes.at(tc.b()), nodes.at(tc.c()));
}

Index FV_Grid::find_N_GHOST_cells()
{
    Index N_GHOST{0};
    for (const auto &tpc : tri_patch_connect_list)
        N_GHOST += tpc.triangles.size();
    return N_GHOST;
}

void FV_Grid::shrink_vectors()
{
    nodes.shrink_to_fit();
    tet_connect.shrink_to_fit();

    for (auto &tpc : tri_patch_connect_list)
        tpc.triangles.shrink_to_fit();

    tri_patch_connect_list.shrink_to_fit();
    cells.cell_volumes.shrink_to_fit();
    cells.centroids.shrink_to_fit();
    faces.cell_indices.shrink_to_fit();
    faces.centroid_to_face_i.shrink_to_fit();
    faces.centroid_to_face_j.shrink_to_fit();
    patches.shrink_to_fit();
}

// void Grid::print_grid(const Config &config) const
// {

//     cout << "CELLS:\n";
//     for (Index i{0}; i < cells.size(); i++)
//     {
//         if (i >= config.get_N_INTERIOR_CELLS())
//             cout << "GHOST, ";
//         cout << i << ": " << cells.at(i) << endl;
//     }

//     cout << "\n\nFACES:\n";
//     for (Index i{0}; i < faces.size(); i++)
//         cout << i << ": " << faces.at(i) << endl;

//     cout << "\n\nPATCHES:\n";
//     for (const auto &patch : patches)
//     {
//         cout << "patch type: " << (int)patch.boundary_type << "\nFIRST FACE: " << patch.FIRST_FACE << "\nN_FACES: " << patch.N_FACES << "\n\n";
//     }
// }

void FV_Grid::print_native_mesh() const
{
    cout << "NODES:\n";
    for (Index i{0}; i < nodes.size(); i++)
    {
        cout << i << ": " << horizontal_string_Vec3(nodes.at(i)) << endl;
    }
    cout << "\n\nTET CONNECTIVITY:\n";
    for (Index i{0}; i < tet_connect.size(); i++)
    {
        TetConnect t = tet_connect.at(i);
        cout << i << ": " << t << endl;
    }
    cout << "\n\nTRIANGLE PATCH CONNECTIVITY:\n";
    for (const auto &tpc : tri_patch_connect_list)
    {
        cout << "BC type: " << tpc.patch_name << endl;
        for (Index i{0}; i < tpc.triangles.size(); i++)
        {
            TriConnect t = tpc.triangles.at(i);
            cout << i << ": " << t << endl;
        }
        cout << "\n\n";
    }
}