#pragma once
#include "../geometry/FV_Grid.hpp"
#include "MetisWrapper.hpp"
#include "Serialization.hpp"

namespace geometry
{
    struct Cellpair
    {
        Index i;
        SignedIndex j;
    };

    class FaceGraph
    {
    protected:
        map<Index, Cellpair> Cellpairs;
        vector<PatchBoundary> patches_bound;

    public:
        const vector<PatchBoundary> &get_patches_ext() const { return patches_bound; }
        void add_patch_bound(PatchBoundary p) { patches_bound.emplace_back(p); }

        bool face_is_in_patch_ext(Index fID, PatchBoundary p) const { return p.FIRST_FACE <= fID && fID < p.FIRST_FACE + p.N_FACES; }
        const map<Index, Cellpair> &get_Cellpairs() const { return Cellpairs; }
        Index size() const { return Cellpairs.size(); }
        Index get_i(Index fID) const
        {
            assert(Cellpairs.count(fID) == 1);
            return Cellpairs.at(fID).i;
        }
        Index get_j(Index fID) const
        {
            assert(Cellpairs.count(fID) == 1);
            return Cellpairs.at(fID).j;
        }
        void add_face(Index i, SignedIndex j) { Cellpairs.emplace(Cellpairs.size(), Cellpair{i, j}); }
        void remove_face(Index fID)
        {
            assert(Cellpairs.count(fID) == 1);
            Cellpairs.erase(fID);
        }
        bool face_is_in_graph(Index fID) const { return Cellpairs.count(fID) == 1; }

        Index find_num_ghost_ext() const
        {
            Index N_GHOST{0};
            for (const auto &p : patches_bound)
                N_GHOST += p.N_FACES;
            assert(N_GHOST > 0);
            return N_GHOST;
        }
    };

    class FaceGraphGlob : public FaceGraph
    {
        const vector<ShortIndex> &part;

    public:
        //     /*The outer vector is a list of all partition patches, and the inner vector is
        //     the global face indices of that patch*/
        //     vector<vector<Index>> contigous_part_faces;

        FaceGraphGlob(const vector<ShortIndex> &part) : part{part} {}

        bool is_ghost_face(Index fID) const
        {
            assert(Cellpairs.count(fID) == 1);
            SignedIndex j = Cellpairs.at(fID).j;
            return j == -1;
        }

        bool is_part_face(Index fID) const
        {
            if (is_ghost_face(fID))
                return false;
            assert(Cellpairs.count(fID) == 1);
            Index i = Cellpairs.at(fID).i;
            SignedIndex j = Cellpairs.at(fID).j;
            return part[i] != part[j];
        }
        bool is_internal_face(Index fID) const
        {
            if (is_ghost_face(fID) || is_part_face(fID))
                return false;
            assert(Cellpairs.count(fID) == 1);
            Index i = Cellpairs.at(fID).i;
            SignedIndex j = Cellpairs.at(fID).j;
            return part[i] == part[j];
        }
    };
    class FaceGraphLoc : public FaceGraph
    {
        Index my_rank;
        vector<PatchInterface> patches_interf; /*Local partition patches seen by a local grid*/

    public:
        FaceGraphLoc(Index my_rank) : my_rank{my_rank} {}
        void add_patch_part(PatchInterface p) { patches_interf.emplace_back(p); }
        const vector<PatchInterface> get_patches_interf() const { return patches_interf; }
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
            for (const auto &p : patches_interf)
                N_GHOST += p.N_FACES;
            assert(N_GHOST > 0);
            return N_GHOST;
        }
    };
    class Utils;
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
        // step 5: create vector<FaceGraphLoc> + vector<vector<Index>> fIDglob2glob + vector<vector<Index>> fIDloc2glob
        // details:
        /*{
            //Need to contain the following: first all internal faces, second all partition faces (grouped), third all
            ghost faces (grouped)
        }*/
        // step 6: create vector<primal_grid_loc> (including surface elements). Input: Primal grid, one of the face graph and utils
        //(Will need nID conversion)
        // Also: vector<Index> eIDglob2loc + vector<Index> eIDloc2glob +
        //                                          vector<vector<Index>> nIDglob2glob + vector<vector<Index>> nIDloc2glob
        // step 7: create vector<FV_grid_loc>. Input: FaceGraphLoc and utils. details:
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
        static vector<FaceGraphLoc> create_local_face_graphs(FaceGraphGlob face_graph,
                                                             Utils &utils);
        /*Step 6*/
        static void create_local_primal_grids(const PrimalGrid &primal_grid,
                                              const vector<FaceGraphLoc> &face_graphs_loc,
                                              const Utils &utils,
                                              vector<unique_ptr<PrimalGrid>> &primal_grids_loc);
        /*Step 7*/
        static void create_local_FV_grids(const vector<FaceGraphLoc> &face_graphs_loc,
                                          vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                          const Utils &utils,
                                          vector<unique_ptr<FV_Grid>> &FV_grids_loc);

        /*Step 8*/
        static void send_recv_grids(Config &config,
                                    vector<unique_ptr<PrimalGrid>> &primal_grids_loc,
                                    vector<unique_ptr<FV_Grid>> &FV_grids_loc,
                                    unique_ptr<PrimalGrid> &primal_grid,
                                    unique_ptr<FV_Grid> &FV_grid);

        /*--------------------------------------------------------------------
        Helper functions
        --------------------------------------------------------------------*/
        /*Reorders faces in a more optimal fashion*/
        static void reorder_face_enitities(Index num_interior_faces,
                                           const vector<PatchInterface> &patches_int,
                                           const vector<PatchBoundary> &patches_ext,
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
        vector<Index> eIDglob2loc_;
        vector<map<Index, Index>> fIDglob2loc_;
        vector<map<Index, Index>> fIDloc2glob_;
        vector<map<Index, Index>> nIDglob2loc_;
        vector<map<Index, Index>> nIDloc2glob_;
        vector<pair<Index, Index>> part2e_range_;
        const vector<ShortIndex> &part;

    public:
        Utils(const vector<ShortIndex> &part) : part{part} {}
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

        const map<Index, Index> &get_nIDglob2loc(ShortIndex r) const { return nIDglob2loc_[r]; }

        void set_eIDglob2loc(vector<Index> &&eIDglob2loc) { eIDglob2loc_ = move(eIDglob2loc); }
        void set_fIDglob2loc(vector<map<Index, Index>> &&fIDglob2loc) { fIDglob2loc_ = move(fIDglob2loc); }
        void set_fIDloc2glob(vector<map<Index, Index>> &&fIDloc2glob) { fIDloc2glob_ = move(fIDloc2glob); }
        void set_nIDglob2loc(vector<map<Index, Index>> &&nIDglob2loc) { nIDglob2loc_ = move(nIDglob2loc); }
        void set_nIDloc2glob(vector<map<Index, Index>> &&nIDloc2glob) { nIDloc2glob_ = move(nIDloc2glob); }
        void set_part2e_range(vector<pair<Index, Index>> &&part2e_range) { part2e_range_ = move(part2e_range); }
    };
}