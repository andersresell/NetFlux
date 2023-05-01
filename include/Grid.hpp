#pragma once

#include "includes.hpp"
#include "Utilities.hpp"
#include "Config.hpp"

namespace geom {

    class Grid{
        //Native mesh
        Vector<Vec3> nodes;
        Vector<TetConnect> tet_connect; 
        Vector<TriPatchConnect> tri_patch_connect;

        //Computational grid
        Vector<Cell> cells;
        Vector<Face> faces;
        Vector<Patch> patches;
        Vector<Vector<Index>> face_indices_from_cell;

    public:    
        
        Grid(const Config& config) {}
        
        void create_grid(Config& config);

        void print_grid(const Config& config) const;
        void print_native_mesh() const;

        const Vector<Vec3>& get_nodes() const {return nodes;}
        const Vector<TetConnect>& get_tet_connect() const {return tet_connect;} 
        const Vector<Cell>& get_cells() const {return cells;}
        const Vector<Face>& get_faces() const {return faces;}

        /*Gets the indices of the faces surrounding cell i*/
        const Vector<Index>& get_surrounding_faces(Index cell_i) const {return face_indices_from_cell[cell_i];}


    private:

        void read_mesh(string mesh_filename);
        void read_c3d_mesh(string mesh_filename);
        void read_su2_mesh(string mesh_filename);
        

        /*------Helper functions for creating grid-------*/

        /*Finds the neighbouring cell j of face ij of cell i. If no neigbour exist (boundary), it returns false*/
        std::pair<Index, bool> find_neigbouring_cell(Index i, 
                                                    TriConnect face_ij, 
                                                    const Vector<TetConnect>& tet_connect) const;
        /*Checks if face ij has been created yet*/
        bool face_ij_created(Index i, Index j) const;

        /*Adds the index of face ij to the correct boundary patch and returns the connectivity of the triangle that ensures normal 
        pointing outwards */
        TriConnect add_face_to_patches(TriConnect t_ij, 
                                       Index ij, 
                                       const Vector<TriPatchConnect>& tri_patch_connect);

        void assign_patch_BC(const Config& config); 

        Tetrahedron tet_from_connect(const TetConnect& tc) const;
        Triangle tri_from_connect(const TriConnect& tc) const;

        /*Adds the indices of neiggbour cell ij to cell i in the face_indices_from_cell structure*/
        void add_face_to_cell_i(Index i, Index ij);

        /*Calling shrink_to_fit on the different members after grid construction to reduce allocated memory*/
        void shrink_vectors();
        
    };

}
