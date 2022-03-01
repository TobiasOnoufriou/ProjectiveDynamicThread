#ifndef _CONVERGE_CUH_
#define _CONVERGE_CUH_

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda.h>
#include <cusparse.h>
//#include "Sparse"

//#include "constraint.h"
//#include <stdio.h>
//#include "constraint.h"

struct CudaConverge {

};

namespace Converge {
	void Converge(double* h_spring, double* h_attach, double* h_pj, double* h_qn1, double* h_b); //This will probably need to be amended for parsing in constraints.
	void CopyMatrixToDevice(const float* matrix, int mRows, int mCols);
}

#endif