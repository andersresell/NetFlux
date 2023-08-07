#include "../../include/parallelization/Communicator.hpp"

using namespace geometry;

InterfaceComm::InterfaceComm(ShortIndex n_vecfields_max,
							 ShortIndex n_gradfields_max,
							 const ShortIndex N_EQS,
							 const geometry::PatchInterface &patch_part,
							 const geometry::Faces &faces)
	: N_EQS{N_EQS}, patch_part{patch_part}, faces{faces}
{
	Index size_vecfield = n_ghost() * N_EQS;
	Index size_gradfield = n_ghost() * N_EQS * N_DIM;
	sendbuf.resize(n_vecfields_max * size_vecfield);
	recvbuf.resize(n_gradfields_max * size_gradfield);
}

void InterfaceComm::send_receive_fields(MPI_Request &send_req, MPI_Request &recv_req)
{

	Index count = sendptr;

	NF_MPI::ISend(sendbuf.data(), count, rank_neigbour(), send_req);

	NF_MPI::IRecv(recvbuf.data(), count, rank_neigbour(), recv_req);
}

PartitionComm::PartitionComm(const ShortIndex n_vecfields_max,
							 const ShortIndex n_gradfields_max,
							 const ShortIndex N_EQS,
							 const FV_Grid &FV_grid)
	: n_vecfields_max{n_vecfields_max}, n_gradfields_max{n_gradfields_max}
{
	const Faces &faces = FV_grid.get_faces();
	const auto &patches_interf = FV_grid.get_patches_interface();
	for (const auto &patch_interf : patches_interf)
		interf_comms.emplace_back(make_unique<InterfaceComm>(n_vecfields_max, n_gradfields_max, N_EQS, patch_interf, faces));
}

void PartitionComm::communicate_ghost_fields()
{

	Vector<MPI_Request> send_requests(num_patches());
	Vector<MPI_Request> recv_requests(num_patches());

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
void PartitionComm::clear()
{
	for (auto &interf_comm : interf_comms)
		interf_comm->clear();
}
