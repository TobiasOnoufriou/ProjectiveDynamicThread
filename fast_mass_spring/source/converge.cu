#include "converge.cuh"


#include <stdio.h>
#include <conio.h>
#include <new>
#include <cmath>

#define SIZE 1024
//ERROR CHECKER CURTOSY OF: https://stackoverflow.com/questions/42180066/cudamemcpy-struct-device-to-host-not-working
#define ERRCHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true, bool wait = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (wait) getch();
		if (abort) exit(code);
	}
}
///// ABOVE CODE CREDITED TO: https://stackoverflow.com/questions/42180066/cudamemcpy-struct-device-to-host-not-working //// 


/// <summary>
/// Functions created below are created by with the aid of the mathematical material.
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <returns></returns>
__device__ double3 operator-(const double3& a, const double3& b) {
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

// localStep is a kernel function that will check each constraint in parrellel then perform convergence to be parsed back to the device.
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


		//Workout p_j's new value 
		//sparseMatrixVectorMultiplication << 1, SIZE >> (cj, pj, d_product, rows); Should be called here but doesn't work.
	}
	__syncthreads();
}

/*sparseMatrixVectorMultiplication should multithred solving matrix multiplication of a sparse matrix and
a vector.*/
__global__ void sparseMatrixVectorMultiplication(double* d_cj, double* p_j, double* d_product, int vRows, int aRow) {
	
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	double product_val = 0.0;
	//m_A(Row) x m_B(column) gives matrix product size matrix.
	//VECTOR rows same size as MATRIX columns
	
		for (int k = 0; k < vRows; k++) {
			product_val += d_cj[(k*vRows)+col] * p_j[k];
		}
		d_product[row*vRows+col] = product_val;
		__syncthreads();
}


// h -> defines host
// d -> defines device
void Converge::Converge(CudaConstraint* h_cj, double* h_b,double* h_pj, double* h_qn1, double* h_pspring, double* h_pattach) {
	double *d_b, *d_pj, *d_qn1, *d_pspring, *d_pattach; //Device memory.
	CudaConstraint* d_cj;
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

/* MatrixMulTest function that handles memory management to the sparseMatrixVectorMultiplication kernel function and will return the product calculated.*/
double* Converge::MatrixMulTest(CudaConstraint* a, double* h_pj, int size) {
	double* d_cj;
	double* d_pj, *d_product;

	double* h_product = (double*)malloc((a->row)*sizeof(double));
	
	int n = a->row * a->col;

	ERRCHECK(cudaMalloc((void**)&d_cj, n * sizeof(double)));
	ERRCHECK(cudaMalloc((void**)&d_product, (a->row) * sizeof(double)));
	ERRCHECK(cudaMalloc((void**)&d_pj, size * sizeof(double)));
	//cudaMalloc((void**)&d_product, (a->row * size) * sizeof(double));
	double* v = a->value;
	ERRCHECK(cudaMemcpy(d_cj,v, n*sizeof(double), cudaMemcpyHostToDevice));
	ERRCHECK(cudaMemcpy(d_pj, h_pj, size * sizeof(double), cudaMemcpyHostToDevice));
	ERRCHECK(cudaMemcpy(d_product, h_product, (a->row) * sizeof(double), cudaMemcpyHostToDevice));
	
	dim3 grid(a->row, a->col);
	sparseMatrixVectorMultiplication<<<grid,1>>>(d_cj, d_pj, d_product, size, a->row);
	
	ERRCHECK(cudaPeekAtLastError());
	ERRCHECK(cudaDeviceSynchronize());

	//cudaMemcpyFromSymbol(&product, d_product, sizeof(product), 0, cudaMemcpyDeviceToHost);
	ERRCHECK(cudaMemcpy(h_product, d_product, (a->row)*sizeof(double), cudaMemcpyDeviceToHost));

	ERRCHECK(cudaFree(d_cj));
	ERRCHECK(cudaFree(d_pj));
	ERRCHECK(cudaFree(d_product));
	//h_product will definetly cause a memory leak.
	return h_product;
}
