#include "converge.cuh"

#define SIZE 1024

/*__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_spring;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_attach;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>* p_j; //May need to be put into the jacobiOnDevice function
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> q_n1;*/



__global__ void localStep(Constraint* m_constraint, VectorX& p_spring, VectorX* p_attach, VectorX* p_j, VectorX* q_n1, VectorX* b) {
	Constraint* cj;
	ScalarType current_strecth; 
	EigenVector3 current_vector;
	int cSize = sizeof m_constraint / sizeof * m_constraint;
	int idx = threadIdx.x;
	int idy = threadIdx.y;
	//Parse in the constraint

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < cSize; i += stride) {
		cj = &m_constraint[i];
		if (cj->constraintType == SPRING) {
			// Work out spring constraint.
			SpringConstraint* sc = (SpringConstraint*) cj;
			current_vector = q_n1->block_vector(sc->GetConstrainedVertexIndex1()) - q_n1->block_vector(sc->GetConstrainedVertexIndex2());
			current_strecth = current_vector.norm() - sc->GetRestLength();
			current_vector = (current_strecth / 2.0) * current_vector.normalized();

			p_j = &p_spring;
			p_j->block_vector(0) = q_n1->block_vector(sc->GetConstrainedVertexIndex1()) - current_vector;
			p_j->block_vector(1) = q_n1->block_vector(sc->GetConstrainedVertexIndex2()) + current_vector;
		}
		if (cj->constraintType == ATTACHMENT) {
			// Work out attachment constraint.
		}
		//cj->m_RHS.applyThisOnTheLeft(*p_j);
		*b += *p_j;
	}

}

// h -> defines host
// d -> defines device
void Converge::Converge(Constraint& cj, VectorX& h_spring, VectorX& h_attach, VectorX& h_pj, VectorX& h_qn1, VectorX& h_b) {

	VectorX* d_b, *d_pj, *d_qn1, *d_pspring, *d_pattach; //Device memory.
	Constraint* d_cj;
	//spring.data 
	//Instead of directly using Eigen use .data and conver it to a float3

	cudaMalloc((void**)&d_cj, sizeof(Constraint));
	cudaMalloc((void**)&d_b, sizeof(VectorX)*h_b.size());
	cudaMalloc((void**)&d_pj, sizeof(VectorX)*h_pj.size());
	//Pspring
	cudaMalloc((void**)&d_pspring, sizeof(VectorX) * h_spring.size());
	//pAttach
	cudaMalloc((void**)&d_pattach, sizeof(VectorX) * h_attach.size());
	//qn1
	cudaMalloc((void**)&d_qn1, sizeof(VectorX) * h_qn1.size());
	

	cudaMemcpy(d_cj, &cj, sizeof(Constraint), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b.data(), sizeof(VectorX)* h_b.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pj, h_pj.data(), sizeof(VectorX) * h_pj.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pspring, h_spring.data(), sizeof(VectorX) * h_spring.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pattach, h_attach.data(), sizeof(VectorX) * h_attach.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_qn1, d_qn1->data(), sizeof(VectorX) * h_qn1.size(), cudaMemcpyHostToDevice);

	localStep<<<1, SIZE >>>(d_cj, 
		d_pspring, 
		d_pattach, 
		d_pj, 
		d_qn1,
		d_b
		);

	cudaDeviceSynchronize();
}
	

//Will need a function that will convert position and velocity to float3
//Returning of p_j will be needed. To be used on the global solver.
