#include "../../include/parallelization/PartComm.hpp"

PartComm::PartComm(const ShortIndex n_vecfields_max,
				   const ShortIndex n_gradfields_max,
				   const ShortIndex N_EQS,
				   const geometry::PatchPart &patch_part,
				   const geometry::Faces &faces)
	: n_vecfields_max{n_vecfields_max}, n_gradfields_max{n_gradfields_max}, N_EQS{N_EQS},
	  patch_part{patch_part}, faces{faces}, size_vecfield{n_ghost() * N_EQS}, size_gradfield{n_ghost() * N_EQS * N_DIM}
{
	sendbuf.resize(n_vecfields_max * size_vecfield);
	recvbuf.resize(n_gradfields_max * size_gradfield);
}

void PartComm::send_and_receive_fields()
{

	assert(sendptr == recvptr);
	MPI_Request send_req;
	MPI_Request recv_req;
	NF_MPI::ISend(sendbuf.data(), sendptr, rank_neigbour(), send_req);

	NF_MPI::IRecv(recvbuf.data(), recvptr, rank_neigbour(), recv_req);
}