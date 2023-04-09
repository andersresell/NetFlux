#pragma once

#include "includes.hpp"
#include "Utilities.hpp"
#include "Config.hpp"



class Grid{
    CellContainer cells;
    vector<Patch> patches;
    FaceContainer faces;




public:    
    
    Grid(const Config& config);
    void create_grid(Config& config);

    void read_mesh(string mesh_filename, 
                   vector<Vec3>& nodes, 
                   vector<Geometry::TetConnectivity>& tet_connect, 
                   vector<Geometry::TriPatchConnectivity>& tri_patch_connect) const;

    /*------Helper functions for creating grid-------*/

    /*Finds the neighbouring cell j of face ij of cell i. If no neigbour exist (boundary), it returns false*/
    std::pair<Index, bool> find_neigbouring_cell(Index i, 
                                                Geometry::TriConnectivity face_ij, 
                                                const vector<Geometry::TetConnectivity>& tet_connect) const;
    /*Checks if face ij has been created yet*/
    bool face_ij_created(Index i, Index j) const;

    /*Adds the index of face ij to the correct boundary patch and returns the connectivity of the triangle that ensures normal 
    pointing outwards */
    Geometry::TriConnectivity add_face_to_patches(Geometry::TriConnectivity t_ij, 
                                                Index ij, 
                                                const vector<Geometry::TriPatchConnectivity>& tri_patch_connect);

};