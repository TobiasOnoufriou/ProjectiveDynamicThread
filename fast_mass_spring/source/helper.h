#pragma once
#include "converge.cuh"


typedef struct {
	float3* m_row;
}VectorX;


inline float3* CreateVectorX(double* matrix, const int row) {
	float3* out = new float3[row];
}



