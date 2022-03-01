#pragma once
class CudaConstraint {
	
	

	

	double* spring_constraints; 
	double* attachement_constraints; //Attachment constraints
	double* m_RHS; //Right hand matrix

};

class CudaSpringConstraint: public CudaConstraint {

	 ///////  Spring Functions  ///////////
	// Getting the Constrained Vertex   // 
	public:
		CudaSpringConstraint(int index1, int index2, double rest_length);
	private: 
		int m_index1, m_index2;
		double m_rest_length;
	public:
		inline unsigned int GetConstrainedVertexIndex1(void) { return m_index1; }
		inline unsigned int GetConstrainedVertexIndex2(void) { return m_index2; }
		inline double GetRestLength(void) { return m_rest_length; }
};

class CudaAttachmentConstraint {

};