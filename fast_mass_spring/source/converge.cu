#include "converge.cuh"

#define SIZE 1024

/*__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_spring;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_attach;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>* p_j; //May need to be put into the jacobiOnDevice function
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> q_n1;*/



__global__ void localStep(double * p_spring, double* p_attach, double* p_j, double* q_n1, double* b ) {
	//Constraint* cj;
	//ScalarType current_strecth; 
	//EigenVector3 current_vector;
	//int cSize = sizeof m_constraint / sizeof * m_constraint;
	int idx = threadIdx.x;
	int idy = threadIdx.y;
	//Parse in the constraint

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	//p_j is just a Matrix.

	//Getting cons

	/*for (int i = index; i < cSize; i += stride) {
		cj = &m_constraint[i];
		if (cj->constraintType == SPRING) {
			// Work out spring constraint.
			SpringConstraint* sc = (SpringConstraint*) cj;
			current_vector = q_n1->block_vector(sc->GetConstrainedVertexIndex1()) - q_n1->block_vector(sc->GetConstrainedVertexIndex2());
			current_strecth = current_vector.norm() - sc->GetRestLength();
			current_vector = (current_strecth / 2.0) * current_vector.normalized();

			p_j = p_spring;
			p_j->block_vector(0) = q_n1->block_vector(sc->GetConstrainedVertexIndex1()) - current_vector;
			p_j->block_vector(1) = q_n1->block_vector(sc->GetConstrainedVertexIndex2()) + current_vector;
		}
		if (cj->constraintType == ATTACHMENT) {
			// Work out attachment constraint.
		}
		//cj->m_RHS.applyThisOnTheLeft(*p_j);
		*b += *p_j;
	}*/

}

// h -> defines host
// d -> defines device
//Return b 
void Converge::Converge(double* h_spring, double* h_attach, double* h_pj, double* h_qn1, double* h_b) {

	double *d_b, *d_pj, *d_qn1, *d_pspring, *d_pattach; //Device memory.
	//Constraint* d_cj;
	//spring.data 
	//Instead of directly using Eigen use .data and conver it to a float3

	d_b = h_b;
	d_pj = h_pj;
	d_qn1 = h_qn1;
	d_pspring = h_spring;
	d_pattach = h_attach;


	//cudaMalloc((void**)&d_cj, sizeof(Constraint));
	cudaMalloc((void**)&d_b, sizeof(double*));
	cudaMalloc((void**)&d_pj, sizeof(double*));
	//Pspring
	cudaMalloc((void**)&d_pspring, sizeof(double*));
	//pAttach
	cudaMalloc((void**)&d_pattach, sizeof(double*));
	//qn1
	cudaMalloc((void**)&d_qn1, sizeof(double*));
	

	//cudaMemcpy(d_cj, &cj, sizeof(Constraint), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, sizeof(double*), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pj, h_pj, sizeof(double*), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pspring, h_spring, sizeof(double*), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pattach, h_attach, sizeof(double*), cudaMemcpyHostToDevice);
	cudaMemcpy(d_qn1, h_qn1, sizeof(double*), cudaMemcpyHostToDevice);

	localStep<<<1, SIZE >>>(
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
