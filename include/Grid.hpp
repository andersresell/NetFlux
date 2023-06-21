#pragma once

#include "includes.hpp"
#include "Utilities.hpp"
#include "Config.hpp"

namespace geom {

    class Grid{
        //Native mesh
        Vector<Vec3> nodes;
        Vector<TetConnect> tet_connect; 
        Vector<Triangle> face_triangles; //Only used for geometry computation, cleared afterwards
        Vector<TriPatchConnect> tri_patch_connect_list;

        //Computational grid
        Vector<Cell> cells;
        Vector<Face> faces;
        Vector<Patch> patches;
        Vector<Vector<Index>> face_indices_from_cell; //Should preferably not be used. More efficient to loop over faces

    public:    
        
        Grid(Config& config);
        

        void print_grid(const Config& config) const;
        void print_native_mesh() const;

        const Vector<Vec3>& get_nodes() const {return nodes;}
        const Vector<TetConnect>& get_tet_connect() const {return tet_connect;} 
        const Vector<Cell>& get_cells() const {return cells;}
        const Vector<Face>& get_faces() const {return faces;}
        const Vector<Patch>& get_patches() const {return patches;}

        /*Gets the indices of the faces surrounding cell i*/
        const Vector<Index>& get_surrounding_faces(Index cell_i) const {return face_indices_from_cell[cell_i];}


    private:


        /*Read mesh file. This populates the:
        - nodes -> nodes of native mesh
        - tet_connect -> tetrahedral elements connectivity
        - tri_patch_connect -> patches of boundary triangles and respective BC types
        */
        void read_mesh(string mesh_filename);
        void read_c3d_mesh(string mesh_filename);
        void read_su2_mesh(string mesh_filename);
        
        /*Creates computational grid from native mesh*/
        void create_grid(Config& config);
        

        /*------Helper functions for creating grid-------*/

        /*Creates interior cells and faces*/
        void create_interior();

        /*loops over the boundary patches and creates faces and associated ghost cells*/
        void create_boundaries(const Config& config);

        //Reorders the (for now interior) faces in an optimal fashion based on the face indices
        void reorder_faces(const Config& config);

        /*Assigns cell centers, boudnary normals, etc*/
        void assign_geometry_properties(const Config& config);

        /*Find the interior cell of the boundary face with a given connectivity*/
        Index find_boundary_face_owner(TriConnect tc);


        /*Finds the neighbouring cell j of face ij of cell i. If no neigbour exist (boundary), it returns false*/
        std::pair<Index, bool> find_neigbouring_cell(Index i, 
                                                    TriConnect face_ij, 
                                                    const Vector<TetConnect>& tet_connect) const;
        /*Checks if face ij has been created yet*/
        bool face_ij_created(Index i, Index j) const;

        Tetrahedron tet_from_connect(const TetConnect& tc) const;
        Triangle tri_from_connect(const TriConnect& tc) const;

        /*Adds the indices of neiggbour cell ij to cell i in the face_indices_from_cell structure*/
        void add_face_to_cell_i(Index i, Index ij);

        /*Used to find the number of ghost cells before this value is set in Config object*/
        Index find_N_GHOST_cells();

        /*Loops through tri_patch_conenct and finds the TriConnect of face ij*/
        //TriConnect find_TriConnect_from_face_index(const Config& config, Index ij)

        /*Calling shrink_to_fit on the different members after grid construction to reduce allocated memory*/
        void shrink_vectors();
        
    };

}
