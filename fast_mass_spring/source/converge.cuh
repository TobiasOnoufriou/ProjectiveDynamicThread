#ifndef _CONVERGE_CUH_
#define _CONVERGE_CUH_

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <thrust/device_vector.h>
#include <device_launch_parameters.h>
#include <cusparse.h>
#include "CudaConstraint.h"


/// <summary>
/// All functions covered in this header file created by Tobias Onoufriou.
/// </summary>
namespace Converge {
	void Converge(CudaConstraint* h_cj, double* h_b, double* h_pj, double* h_qn1, double* h_pspring, double* h_pattach); //This will probably need to be amended for parsing in constraints.
	double* MatrixMulTest(CudaConstraint* a, double* p_j, int size);
	//TODO: multiplication on left hand side of matrix like applyThisOnTheLeft.
}
#endif
