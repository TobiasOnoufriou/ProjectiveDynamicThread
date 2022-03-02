#include "CudaConstraint.h"

CudaAttachmentConstraint::CudaAttachmentConstraint(double* fixedPoint, unsigned int vertex_index) {
	this->m_fixed_point = fixedPoint;
	this->m_vertex_index = vertex_index;
}


CudaSpringConstraint::CudaSpringConstraint(unsigned int index1, unsigned int index2, double rest_length) {
	this->m_index1 = index1;
	this->m_index2 = index2;
	this->m_rest_length = rest_length;
}