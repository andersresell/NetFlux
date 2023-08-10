#include "../../include/parallelization/Communicator.hpp"

using namespace geometry;

InterfaceComm::InterfaceComm(const geometry::PatchInterface &patch_part,
							 const geometry::Faces &faces)
	: patch_part{patch_part}, faces{faces} {}
//
// 	Index size_vecfield = n_ghost() * N_EQS;
// 	Index size_gradfield = n_ghost() * N_EQS * N_DIM;
// 	sendbuf.resize(n_vecfields_max * size_vecfield);
// 	recvbuf.resize(n_gradfields_max * size_gradfield);

void InterfaceComm::send_receive_fields(MPI_Request &send_req, MPI_Request &recv_req)
{

	Index count = sendptr;

	NF_MPI::ISend(sendbuf.data(), count, rank_neigbour(), send_req);

	NF_MPI::IRecv(recvbuf.data(), count, rank_neigbour(), recv_req);
}

PartitionComm::PartitionComm(const Config &config, const geometry::FV_Grid &FV_grid)
	: N_CELLS_TOT{config.get_N_CELLS_TOT()}
{
	assert(FV_grid.get_cells().size() == config.get_N_CELLS_TOT());
	const Faces &faces = FV_grid.get_faces();
	const auto &patches_interf = FV_grid.get_patches_interface();
	for (const auto &patch_interf : patches_interf)
		interf_comms.emplace_back(make_unique<InterfaceComm>(patch_interf, faces));
}

void PartitionComm::communicate_ghost_fields()
{

	vector<MPI_Request> send_requests(num_patches());
	vector<MPI_Request> recv_requests(num_patches());

	for (ShortIndex i{0}; i < interf_comms.size(); i++)
	{
		interf_comms[i]->send_receive_fields(send_requests[i], recv_requests[i]);
	}
	for (ShortIndex i{0}; i < interf_comms.size(); i++)
	{
		NF_MPI::Wait(send_requests[i]);
		NF_MPI::Wait(recv_requests[i]);
	}
}

void PartitionComm::communicate_interface_ghost_centroids(Vector<Vec3> &centroids)
{
	if (get_max_size_cell() < N_DIM)
		set_max_size_cell(N_DIM);

	clear();
	pack_Vec3_field(centroids);
	communicate_ghost_fields();
	unpack_Vec3_field(centroids);
}
void PartitionComm::clear()
{
	for (auto &interf_comm : interf_comms)
		interf_comm->clear();
}
