#ifndef _CONVERGE_CUH_
#define _CONVERGE_CUH_

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda.h>
#include <cusparse.h>
#include "CudaConstraint.h"

namespace Converge {
	void Converge(CudaConstraint* c, double* p_j, double* q_n1, double* p_spring, double* p_attach); //This will probably need to be amended for parsing in constraints.
}

#endif