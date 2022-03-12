#include "converge.cuh"


#define SIZE 1024

/*__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_spring;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> p_attach;
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>* p_j; //May need to be put into the jacobiOnDevice function
__device__ Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> q_n1;*/


__device__ double* d_product;

/*__device__ double3 operator-(const double3& a, const double3& b) {
	return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ double3 operator+(const double3& a, const double3& b) {
	return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__device__ double3 operator*(const double& a, const double3& b) {
	return make_double3(a * b.x, a * b.y, a * b.z);
}

__device__ double norm(double3 a) {
	return sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
}

__device__ double3 normalise(double3 a){
	double3 out;
	float s = sqrt((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
	out = make_double3(a.x / s, a.y / s, a.z / s);
	return out;
}

__device__ double* double3ToDouble(double3 a) {
	double* out = (double*)malloc(3);
	out[0] = a.x;
	out[1] = a.y;
	out[2] = a.z;
	return out;
}

__device__ void ajoinVectorArray(double* a, double* b, double* c) {
	for (int i = 0; i < 6; i++) {
		if (i < 3) {
			c[i] = a[i];
		}
		else {
			c[i] = b[i];
		}
	}
}

__global__ void localStep(CudaConstraint* d_cj, double* p_spring, double* p_attach, double* q_n1, double* b ) {
	CudaConstraint* cj;
	CudaSpringConstraint* scj;
	CudaAttachmentConstraint* acj;
	int rows;
	double3* pj;
	
	double current_strecth; 
	double3 current_vector;

	//Parse in the constraint
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	
	__syncthreads();
	
	//p_j is just a Matrix.
	if (index < sizeof(d_cj)) {
		cj = &d_cj[index];

		if (cj->constraint == SPR) {
			scj = (CudaSpringConstraint*)cj;
			double3 V1 = make_double3(q_n1[scj->m_index1 * 1], q_n1[scj->m_index1 * 2], q_n1[scj->m_index1 * 3]);
			double3 V2 = make_double3(q_n1[scj->m_index2 * 1], q_n1[scj->m_index2 * 2], q_n1[scj->m_index2 * 3]);
			current_vector = V1 - V2;
			current_strecth = norm(current_vector) - scj->m_rest_length;
			current_vector = (current_strecth / 2.0) * normalise(current_vector);
			
			double3 tPj[2];
			tPj[0] = V1 - current_vector;
			tPj[1] = V2 + current_vector;
			
			pj = tPj;
			rows = 6;
			//Call "applyToLeft function here
		}

		else if (cj->constraint == ATT) {
			acj = (CudaAttachmentConstraint*)cj;

			double3 tPj[1];
			tPj[0] = make_double3(acj->m_fixed_point[0], acj->m_fixed_point[1], acj->m_fixed_point[2]);
			rows = 3;
			pj = tPj;
		}
		//TODO: Add sparse matrix applying from left. Which seems just to multiply on the right e.g m_RHS * p_j
		//TODO: Matrix multiplication with b to be returned and sovled.



		double* d_product = (double*)malloc(nRows * d_cj.col);
		//Workout p_j's new value 
		//sparseMatrixVectorMultiplication << 1, SIZE >> (cj, pj, d_product, rows);

	}
	
	__syncthreads();
}*/

//Attempt at adding this first before adding it to the bigger system.
/*sparseMatrixVectorMultiplication should multithred solving matrix multiplication of a sparse matrix and
a vector.*/
__global__ void sparseMatrixVectorMultiplication(CudaConstraint* d_cj, double* p_j, int nRows) {
	int i = blockDim.y * blockIdx.y + threadIdx.y;
	int j = blockDim.x * blockIdx.x + threadIdx.x;
	double product_val = 0;
	//m_A(Row) x m_B(column) gives matrix product size
	//VECTOR rows same size as MATRIX columns
	for (int k = 0; k < nRows; k++){
		product_val += d_cj->value[k * d_cj->row + j] * p_j[i * d_cj->row + k];
	}
	d_product[i * nRows + j] = product_val;
}


// h -> defines host
// d -> defines device
void Converge::Converge(CudaConstraint* h_cj, double* h_b,double* h_pj, double* h_qn1, double* h_pspring, double* h_pattach) {

	double *d_b, *d_pj, *d_qn1, *d_pspring, *d_pattach; //Device memory.
	CudaConstraint* d_cj;
	//Instead of directly using Eigen use .data and conver it to a float3
	// Conversion Dense to Sparse using cusparseDenseToSparse


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

/*MatrixMulTest function that handles memory management to the sparseMatrixVectorMultiplication kernel function and will return the product calculated.*/
double* Converge::MatrixMulTest(CudaConstraint* a, double* h_pj, double size) {
	CudaConstraint* d_cj;
	double* d_pj, *d_product;

	double* h_product;// = (double*)malloc(a->row * size);
	dim3 grid(1, 1), block(a->row, 1);

	cudaMalloc((void**)&d_cj, sizeof(CudaConstraint*));

	cudaMalloc((void**)&d_pj, size * sizeof(double));
	//cudaMalloc((void**)&d_product, (a->row * size) * sizeof(double));

	cudaMemcpy(&a, d_cj, sizeof(CudaConstraint*), cudaMemcpyHostToDevice);
	cudaMemcpy(&h_pj, d_pj, size * sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(&h_product, d_product, (a->row * size) * sizeof(double), cudaMemcpyHostToDevice);

	sparseMatrixVectorMultiplication<<<grid, block>>>(d_cj, d_pj, size);
	cudaDeviceSynchronize();

	//cudaMemcpyFromSymbol(&product, d_product, sizeof(product), 0, cudaMemcpyDeviceToHost);
	cudaMemcpy(&d_product, h_product, (a->row * size)*sizeof(double), cudaMemcpyDeviceToHost);
	//std::cout << sizeof(h_product)/sizeof(h_product[0]) << std::endl;

	cudaFree(d_cj);
	cudaFree(d_pj);
	cudaFree(d_product);
	
	return h_product;
}
	
void Converge::ConvertDenseToCuSparse(CudaConstraint* h_cj) {
}
//Will need a function that will convert position and velocity to float3
//Returning of p_j will be needed. To be used on the global solver.
