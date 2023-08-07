#pragma once
#include "../geometry/FV_Grid.hpp"
#include "Metis_Wrapper.hpp"
#include "Serialization.hpp"

namespace geometry
{
    struct Cellpair
    {
        Index i;
        long int j;
    };

    class FaceGraph
    {
    protected:
        map<Index, Cellpair> cellpairs;
        Vector<PatchBoundary> patches_ext;

    public:
        const Vector<PatchBoundary> &get_patches_ext() const { return patches_ext; }
        void add_patch_ext(PatchBoundary p) { patches_ext.emplace_back(p); }

        bool face_is_in_patch_ext(Index fID, PatchBoundary p) const { return p.FIRST_FACE <= fID && fID < p.FIRST_FACE + p.N_FACES; }
        const map<Index, Cellpair> &get_cellpairs() const { return cellpairs; }
        Index size() const { return cellpairs.size(); }
        Index get_i(Index fID) const
        {
            assert(cellpairs.count(fID) == 1);
            return cellpairs.at(fID).i;
        }
        Index get_j(Index fID) const
        {
            assert(cellpairs.count(fID) == 1);
            return cellpairs.at(fID).j;
        }
        void add_face(Index i, Index j) { cellpairs.emplace(cellpairs.size(), Cellpair{i, j}); }
        void remove_face(Index fID)
        {
            assert(cellpairs.count(fID) == 1);
            cellpairs.erase(fID);
        }
        bool face_is_in_graph(Index fID) const { return cellpairs.count(fID) == 1; }

        Index find_num_ghost_ext() const
        {
            Index N_GHOST{0};
            for (const auto &p : patches_ext)
                N_GHOST += p.N_FACES;
            assert(N_GHOST > 0);
            return N_GHOST;
        }
    };

    class FaceGraphGlob : public FaceGraph
    {
        const Vector<ShortIndex> &part;

    public:
        //     /*The outer vector is a list of all partition patches, and the inner vector is
        //     the global face indices of that patch*/
        //     Vector<Vector<Index>> contigous_part_faces;

        FaceGraphGlob(const Vector<ShortIndex> &part) : part{part} {}

        bool is_ghost_face(Index fID) const
        {
            assert(cellpairs.count(fID) == 1);
            Index i = cellpairs.at(fID).i;
            Index j = cellpairs.at(fID).j;
            return i == -1 || j == -1;
        }

        bool is_part_face(Index fID) const
        {
            if (is_ghost_face(fID))
                return false;
            assert(cellpairs.count(fID) == 1);
            Index i = cellpairs.at(fID).i;
            Index j = cellpairs.at(fID).j;
            return part[i] != part[j];
        }
        bool is_internal_face(Index fID) const
        {
            if (is_ghost_face(fID) || is_part_face(fID))
                return false;
            assert(cellpairs.count(fID) == 1);
            Index i = cellpairs.at(fID).i;
            Index j = cellpairs.at(fID).j;
            return part[i] == part[j];
        }
    };
    class FaceGraphLoc : public FaceGraph
    {
        Index my_rank;
        Vector<PatchInterface> patches_int; /*Local partition patches seen by a local grid*/

    public:
        FaceGraphLoc(Index my_rank) : my_rank{my_rank} {}
        void add_patch_part(PatchInterface p) { patches_int.emplace_back(p); }
        const Vector<PatchInterface> get_patches_int() const { return patches_int; }
        //  void add_PatchBoundary(BoundaryType bt, Index begin)
        // {
        //     PatchBoundary p;
        //     p.boundary_type = bt;
        //     p.FIRST_FACE = begin;
        //     patches.push_back(p);
        // }
        // void set_PatchBoundary_end(Index end)
        // {
        //     patches.back().N_FACES = end - patches.back().FIRST_FACE;
        // }

        Index find_num_ghost_part() const
        {
            Index N_GHOST{0};
            for (const auto &p : patches_int)
                N_GHOST += p.N_FACES;
            assert(N_GHOST > 0);
            return N_GHOST;
        }
    };

    class GridCreator
    {
    public:
        static void create_partitioned_grids(Config &config,
                                             unique_ptr<PrimalGrid> &primal_grid_glob,
                                             unique_ptr<PrimalGrid> &primal_grid,
                                             unique_ptr<FV_Grid> &FV_grid);

        // step 1: read native mesh and create primal_grid_glob
        // step 2: run metis and create part. enables eidloc2glob
        // step 3: reorder primal_grid_glob, grouping each rank's elements together. sets eIDglob2loc
        // step 4: create FaceGraphGlob + FaceElementsGlob
        // details:
        /*{
            Need to include details about ghost cells and boundaries.
            Should contain
        }*/
        // step 5: create Vector<FaceGraphLoc> + vector<vector<Index>> fIDglob2glob + vector<vector<Index>> fIDloc2glob
        // details:
        /*{
            //Need to contain the following: first all internal faces, second all partition faces (grouped), third all
            ghost faces (grouped)
        }*/
        // step 6: create Vector<primal_grid_loc> (including surface elements). Input: Primal grid, one of the face graph and utils
        //(Will need nID conversion)
        // Also: vector<Index> eIDglob2loc + vector<Index> eIDloc2glob +
        //                                          vector<vector<Index>> nIDglob2glob + vector<vector<Index>> nIDloc2glob
        // step 7: create Vector<FV_grid_loc>. Input: FaceGraphLoc and utils. details:
        /*{
            //Create each one simply from FaceGraphLoc, get the correct faces by using e_faces, fIDloc2glob and nIDglob2loc
        }*/

        // Step 8: Send and receive all grids from rank 0 to the other procs

        /*Question for later: What is the requirements to the local indices of the PatchBoundary faces, given that packing is used
        (each PatchBoundary of each rank has a sendbuf and recvbuf). Conclusion is no requirements, but */
        /*Comment: I realize that I haven't accounted for the possibility that two ranks can share more than one interface.
        I will ignore this for now*/

    private:
        /*Step 3*/
        static void reorder_global_grid(PrimalGrid &primal_grid,
                                        Utils &utils);

        /*Step 4*/
        static void create_global_face_graph_and_face_elements(const Config &config,
                                                               PrimalGrid &primal_grid,
                                                               FaceGraphGlob &face_graph,
                                                               Utils &utils);
        /*Step 5*/
        static Vector<FaceGraphLoc> create_local_face_graphs(FaceGraphGlob face_graph,
                                                             Utils &utils);
        /*Step 6*/
        static void create_local_primal_grids(const PrimalGrid &primal_grid,
                                              const Vector<FaceGraphLoc> &face_graphs_loc,
                                              const Utils &utils,
                                              Vector<unique_ptr<PrimalGrid>> &primal_grids_loc);
        /*Step 7*/
        static void create_local_FV_grids(const Vector<FaceGraphLoc> &face_graphs_loc,
                                          Vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                          const Utils &utils,
                                          Vector<unique_ptr<FV_Grid>> &FV_grids_loc);

        /*Step 8*/
        static void send_recv_grids(Config &config,
                                    Vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                    Vector<unique_ptr<FV_Grid>> &FV_grids_loc,
                                    unique_ptr<PrimalGrid> &primal_grid,
                                    unique_ptr<FV_Grid> &FV_grid);

        /*--------------------------------------------------------------------
        Helper functions
        --------------------------------------------------------------------*/
        /*Reorders faces in a more optimal fashion*/
        static void reorder_face_enitities(Index num_interior_faces,
                                           const Vector<PatchInterface> &patches_int,
                                           const Vector<PatchBoundary> &patches_ext,
                                           Faces &faces,
                                           Elements &face_elements);

        static void set_global_config_data(Config &config,
                                           unique_ptr<PrimalGrid> &primal_grid_glob,
                                           const FaceGraphGlob &face_graph_glob);

        static void set_config_grid_data_local(Config &config,
                                               unique_ptr<PrimalGrid> &primal_grid,
                                               unique_ptr<FV_Grid> &FV_grid);
    };

    class Utils
    {
        Vector<Index> eIDglob2loc_;
        Vector<map<Index, Index>> fIDglob2loc_;
        Vector<map<Index, Index>> fIDloc2glob_;
        Vector<map<Index, Index>> nIDglob2loc_;
        Vector<map<Index, Index>> nIDloc2glob_;
        Vector<pair<Index, Index>> part2e_range_;
        const Vector<ShortIndex> &part;

    public:
        Utils(const Vector<ShortIndex> &part) : part{part} {}
        Index eIDglob2loc(Index eIDglob) const
        {
            assert(eIDglob2loc_.size() > 0);
            return eIDglob2loc_[eIDglob];
        };
        Index eIDloc2glob(ShortIndex r, Index eIDloc) const
        {
            assert(part2e_range_.size() == NF_MPI::get_size());
            return eIDloc + part2e_range_[r].first;
        };
        Index fIDglob2loc(ShortIndex r, Index fIDglob) const
        {
            assert(fIDglob2loc_.size() == NF_MPI::get_size());
            assert(fIDglob2loc_[r].count(fIDglob) == 1);
            return fIDglob2loc_[r].at(fIDglob);
        };
        Index fIDloc2glob(ShortIndex r, Index fIDloc) const
        {
            assert(fIDloc2glob_.size() == NF_MPI::get_size());
            assert(fIDloc2glob_[r].count(fIDloc) == 1);
            return fIDloc2glob_[r].at(fIDloc);
        };
        Index nIDglob2loc(ShortIndex r, Index nIDglob) const
        {
            assert(nIDglob2loc_.size() == NF_MPI::get_size());
            assert(nIDglob2loc_[r].count(nIDglob) == 1);
            return nIDglob2loc_[r].at(nIDglob);
        }
        Index nIDloc2glob(ShortIndex r, Index nIDloc) const
        {
            assert(nIDloc2glob_.size() == NF_MPI::get_size());
            assert(nIDloc2glob_[r].count(nIDloc) == 1);
            return nIDloc2glob_[r].at(nIDloc);
        }
        Index part2e_range_begin(ShortIndex r) const
        {
            assert(nIDloc2glob_.size() == NF_MPI::get_size());
            return part2e_range_[r].first;
        }

        Index part2e_range_end(ShortIndex r) const
        {
            assert(nIDloc2glob_.size() == NF_MPI::get_size());
            return part2e_range_[r].second;
        }
        ShortIndex e2r(Index eID) const { return part[eID]; }

        Vector<pair<Index, Index>> part2e_range_;

        const map<Index, Index> &get_nIDglob2loc(ShortIndex r) const { return nIDglob2loc_[r]; }

        void set_eIDglob2loc(Vector<Index> &&eIDglob2loc) { eIDglob2loc_ = move(eIDglob2loc); }
        void set_fIDglob2loc(Vector<map<Index, Index>> &&fIDglob2loc) { fIDglob2loc_ = move(fIDglob2loc); }
        void set_fIDloc2glob(Vector<map<Index, Index>> &&fIDloc2glob) { fIDloc2glob_ = move(fIDloc2glob); }
        void set_nIDglob2loc(Vector<map<Index, Index>> &&nIDglob2loc) { nIDglob2loc_ = move(nIDglob2loc); }
        void set_nIDloc2glob(Vector<map<Index, Index>> &&nIDloc2glob) { nIDloc2glob_ = move(nIDloc2glob); }
        void set_part2e_range(Vector<pair<Index, Index>> &&part2e_range) { part2e_range_ = move(part2e_range); }
    };
}