#ifndef _CONVERGE_CUH_
#define _CONVERGE_CUH_

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda.h>

#include "constraint.h"
#include <stdio.h>
//#include "constraint.h"

namespace Converge {
	void Converge(Constraint& cj, VectorX& h_spring, VectorX& h_attach, VectorX& h_pj, VectorX& h_qn1, VectorX& h_b); //This will probably need to be amended for parsing in constraints.
	void CopyMatrixToDevice(const float* matrix, int mRows, int mCols);
}

#endif