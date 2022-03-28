#pragma once

// Attachment and Spring constraint
enum {ATT, SPR};

class CudaConstraint {
public:
	//Sparse Matrix
	int num_non0;
	int num_outer;

	int row, col;
	double* value; //Matrix

	int constraint;
};

class CudaSpringConstraint: public CudaConstraint {

	 ///////  Spring Functions  ///////////
	// Getting the Constrained Vertex   // 
	public:
		CudaSpringConstraint(unsigned int index1, unsigned int index2, double rest_length);
	public: 
		unsigned int m_index1, m_index2;
		double m_rest_length;
	public:
		inline unsigned int GetConstrainedVertexIndex1(void) { return m_index1; }
		inline unsigned int GetConstrainedVertexIndex2(void) { return m_index2; }
		inline double GetRestLength(void) { return m_rest_length; }
};

class CudaAttachmentConstraint: public CudaConstraint {
public:
	CudaAttachmentConstraint(double* fixedPoint, unsigned int vertex_index);
public:
	double* m_fixed_point;
	unsigned int m_vertex_index;
public:
	inline double* GetAttachmentFixedPoint() { return m_fixed_point; }
};