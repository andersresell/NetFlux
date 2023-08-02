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
        Vector<Patch> patches;

    public:
        const Vector<Patch> &get_patches() const { return patches; }
        void add_patch(Patch p) { patches.emplace_back(p); }

        bool face_is_in_patch(Index fID, Patch p) const { return p.FIRST_FACE <= fID && fID < p.FIRST_FACE + p.N_FACES; }
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
    };

    /*--------------------------------------------------------------------
    Contains:
    All global faces with global cell pairs
    All external patches
    Ordering: internal patches first, then each external patch.
    --------------------------------------------------------------------*/
    class FaceGraphGlob : public FaceGraph
    {
        map<Index, Index> part_faces; /*Pointers to the face ID's that separates two partitions*/

        const Vector<ShortIndex> &part;

    public:
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
    /*--------------------------------------------------------------------
    Local face graoh for rank r
    Contains:
    * All local faces with local cell indices
    * A vector of partitions for all ranks (size=num_procs),
      the ones not connecting or its own rank will be empty
    * All boundary patches of its own rank.
    --------------------------------------------------------------------*/
    class FaceGraphLoc : public FaceGraph
    {
        Index my_rank;
        Vector<pair<Index, Index>> part_begin_end;
        Vector<Patch> patches;

    public:
        FaceGraphLoc(Index my_rank) : my_rank{my_rank}, part_begin_end{NF_MPI::get_size()} {}

        void set_part_begin(ShortIndex r, Index begin)
        {
            assert(r < NF_MPI::get_size());
            part_begin_end[r].first = begin;
        }
        void set_part_end(ShortIndex r, Index end)
        {
            assert(r < NF_MPI::get_size());
            part_begin_end[r].second = end;
        }
        void add_patch(BoundaryType bt, Index begin)
        {
            Patch p;
            p.boundary_type = bt;
            p.FIRST_FACE = begin;
            patches.push_back(p);
        }

        void set_patch_end(Index end)
        {
            patches.back().N_FACES = end - patches.back().FIRST_FACE;
        }
    };

    class GridCreatorNew
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

        /*Question for later: What is the requirements to the local indices of the patch faces, given that packing is used
        (each patch of each rank has a sendbuf and recvbuf). Conclusion is no requirements, but */
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
        static Vector<PrimalGrid> create_local_primal_grids(const PrimalGrid &primal_grid,
                                                            const Vector<FaceGraphLoc> &face_graphs_loc,
                                                            Utils &utils);
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