#include "converge.cuh"

#define SIZE 1024

/*__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_spring;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_attach;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>* p_j; //May need to be put into the jacobiOnDevice function
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> q_n1;*/


__device__ float3 operator-(const float3& a, const float3& b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ float3 operator+(const float3& a, const float3& b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ float3 operator*(const double& a, const float3& b) {
	return make_float3(a * b.x, a * b.y, a * b.z);
}

__device__ double norm(float3 a) {
	return sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
}

__device__ float3 normalise(float3 a){
	float3 out;
	float s = sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
	out = make_float3(a.x / s, a.y / s, a.z / s);
	return out;
}

__global__ void localStep(CudaConstraint* d_cj, double* p_spring, double* p_attach, double* q_n1, double* b ) {
	CudaConstraint* cj;
	CudaSpringConstraint* scj;
	CudaAttachmentConstraint* acj;
	
	float3* pj;
	
	double current_strecth; 
	float3 current_vector;

	//Parse in the constraint
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	
	__syncthreads();
	
	//p_j is just a Matrix.

	if (index < sizeof(d_cj)) {
		cj = &d_cj[index];

		if (cj->constraint == SPR) {
			scj = (CudaSpringConstraint*)cj;
			float3 V1 = make_float3(q_n1[scj->m_index1 * 1], q_n1[scj->m_index1 * 2], q_n1[scj->m_index1 * 3]);
			float3 V2 = make_float3(q_n1[scj->m_index2 * 1], q_n1[scj->m_index2 * 2], q_n1[scj->m_index2 * 3]);
			current_vector = V1 - V2;
			current_strecth = norm(current_vector) - scj->m_rest_length;
			current_vector = (current_strecth / 2.0) * normalise(current_vector);
			
			float3 tPj[2];
			tPj[0] = V1 - current_vector;
			tPj[1] = V2 - current_vector;

			pj = tPj;
		}

		else if (cj->constraint == ATT) {
			acj = (CudaAttachmentConstraint*)cj;

			float3 tPj;
			tPj = make_float3(acj->m_fixed_point[0], acj->m_fixed_point[1], acj->m_fixed_point[2]);
		}

		//TODO: Add sparse matrix applying from left. Which seems just to multiply on the right e.g m_RHS * p_j
		//TODO: Matrix multiplication with b to be returned and sovled.
	}
	
	__syncthreads();
}

// h -> defines host
// d -> defines device
void Converge::Converge(CudaConstraint* h_cj, double* h_b,double* h_pj, double* h_qn1, double* h_pspring, double* h_pattach) {

	double *d_b, *d_pj, *d_qn1, *d_pspring, *d_pattach; //Device memory.
	CudaConstraint* d_cj;
	//Instead of directly using Eigen use .data and conver it to a float3

	cusparseMatDescr_t Adescr = 0;

	

	//d_b = h_b;
	//d_pj = h_pj;
	//d_qn1 = h_qn1;
	//d_pspring = h_spring;
	//d_pattach = h_attach;

	//cudaMalloc((void**)&d_cj, sizeof(Constraint));
	cudaMalloc((void**)&d_b, sizeof(double*));
	cudaMalloc((void**)&d_pj, sizeof(double*));
	//Pspring
	cudaMalloc((void**)&d_pspring, sizeof(double*));
	//pAttach
	cudaMalloc((void**)&d_pattach, sizeof(double*));
	//qn1
	cudaMalloc((void**)&d_qn1, sizeof(double*));
	

	cudaMemcpy(d_cj, &h_cj, sizeof(CudaConstraint), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, &h_b, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pj, &h_pj, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pspring, &h_pspring, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_pattach, &h_pattach, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_qn1, &h_qn1, sizeof(double), cudaMemcpyHostToDevice);

	/*localStep<<<1, SIZE >>>(
		d_pspring,
		d_pattach, 
		d_pj, 
		d_qn1,
		d_b
		);*/
	
	cudaDeviceSynchronize();
}
	

//Will need a function that will convert position and velocity to float3
//Returning of p_j will be needed. To be used on the global solver.
