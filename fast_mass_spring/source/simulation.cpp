// ---------------------------------------------------------------------------------//
// Copyright (c) 2013, Regents of the University of Pennsylvania                    //
// All rights reserved.                                                             //
//                                                                                  //
// Redistribution and use in source and binary forms, with or without               //
// modification, are permitted provided that the following conditions are met:      //
//     * Redistributions of source code must retain the above copyright             //
//       notice, this list of conditions and the following disclaimer.              //
//     * Redistributions in binary form must reproduce the above copyright          //
//       notice, this list of conditions and the following disclaimer in the        //
//       documentation and/or other materials provided with the distribution.       //
//     * Neither the name of the <organization> nor the                             //
//       names of its contributors may be used to endorse or promote products       //
//       derived from this software without specific prior written permission.      //
//                                                                                  //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  //
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    //
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           //
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY               //
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       //
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     //
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      //
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       //
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     //
//                                                                                  //
// Contact Tiantian Liu (ltt1598@gmail.com) if you have any questions.              //
//----------------------------------------------------------------------------------//

#pragma warning( disable : 4996)
#include <omp.h>
#include <exception>

#include "simulation.h"
#include "timer_wrapper.h"

#include "simpleTimer.h"

TimerWrapper g_integration_timer;

Simulation::Simulation()
{
	////////////////////////////////////////////////////
	// setup A matrix for attachment constraint
	////////////////////////////////////////////////////

	std::vector<SparseMatrixTriplet> att_triplets;

	float av1 = 1;

	att_triplets.push_back( SparseMatrixTriplet( 0, 0, av1 ) );
	att_triplets.push_back( SparseMatrixTriplet( 1, 1, av1 ) );
	att_triplets.push_back( SparseMatrixTriplet( 2, 2, av1 ) );

	m_A_attachment.resize(3,3);
	m_A_attachment.setFromTriplets( att_triplets.begin(), att_triplets.end() );


	////////////////////////////////////////////////////
	// setup A matrix for spring constraint
	////////////////////////////////////////////////////

	std::vector<SparseMatrixTriplet> spr_triplets;

	float sv1 = 0.5;
	float sv2 = -0.5;

	spr_triplets.push_back( SparseMatrixTriplet( 0, 0, sv1 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 1, 1, sv1 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 2, 2, sv1 ) );
	
	spr_triplets.push_back( SparseMatrixTriplet( 0, 3, sv2 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 1, 4, sv2 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 2, 5, sv2 ) );
	
	spr_triplets.push_back( SparseMatrixTriplet( 3, 0, sv2 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 4, 1, sv2 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 5, 2, sv2 ) );

	spr_triplets.push_back( SparseMatrixTriplet( 3, 3, sv1 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 4, 4, sv1 ) );
	spr_triplets.push_back( SparseMatrixTriplet( 5, 5, sv1 ) );

	m_A_spring.resize(6,6);
	m_A_spring.setFromTriplets( spr_triplets.begin(), spr_triplets.end() );


	////////////////////////////////////////////////////
	// setup A matrix for tet constraint
	////////////////////////////////////////////////////

	std::vector<SparseMatrixTriplet> tet_triplets;

	float v1 =  (2.0f / 3.0f);
	float v2 = -(1.0f / 3.0f);

	// block 1,1 = v1 * I4									// block 1,2 = v2 * I4									// block 1,3 = v2 * I4
	tet_triplets.push_back(SparseMatrixTriplet( 0, 0, v1));	tet_triplets.push_back(SparseMatrixTriplet( 0, 4, v2));	tet_triplets.push_back(SparseMatrixTriplet( 0,  8, v2));
	tet_triplets.push_back(SparseMatrixTriplet( 1, 1, v1));	tet_triplets.push_back(SparseMatrixTriplet( 1, 5, v2));	tet_triplets.push_back(SparseMatrixTriplet( 1,  9, v2));
	tet_triplets.push_back(SparseMatrixTriplet( 2, 2, v1));	tet_triplets.push_back(SparseMatrixTriplet( 2, 6, v2));	tet_triplets.push_back(SparseMatrixTriplet( 2, 10, v2));
	tet_triplets.push_back(SparseMatrixTriplet( 3, 3, v1));	tet_triplets.push_back(SparseMatrixTriplet( 3, 7, v2));	tet_triplets.push_back(SparseMatrixTriplet( 3, 11, v2));

	// block 2,1 = v2 * I4									// block 2,2 = v1 * I4									// block 2,3 = v2 * I4
	tet_triplets.push_back(SparseMatrixTriplet( 4, 0, v2));	tet_triplets.push_back(SparseMatrixTriplet( 4, 4, v1));	tet_triplets.push_back(SparseMatrixTriplet( 4,  8, v2));
	tet_triplets.push_back(SparseMatrixTriplet( 5, 1, v2));	tet_triplets.push_back(SparseMatrixTriplet( 5, 5, v1));	tet_triplets.push_back(SparseMatrixTriplet( 5,  9, v2));
	tet_triplets.push_back(SparseMatrixTriplet( 6, 2, v2));	tet_triplets.push_back(SparseMatrixTriplet( 6, 6, v1));	tet_triplets.push_back(SparseMatrixTriplet( 6, 10, v2));
	tet_triplets.push_back(SparseMatrixTriplet( 7, 3, v2));	tet_triplets.push_back(SparseMatrixTriplet( 7, 7, v1));	tet_triplets.push_back(SparseMatrixTriplet( 7, 11, v2));
	
	// block 3,1 = v2 * I4									// block 3,2 = v2 * I4									// block 3,3 = v1 * I4
	tet_triplets.push_back(SparseMatrixTriplet( 8, 0, v2));	tet_triplets.push_back(SparseMatrixTriplet( 8, 4, v2));	tet_triplets.push_back(SparseMatrixTriplet( 8,  8, v1));
	tet_triplets.push_back(SparseMatrixTriplet( 9, 1, v2));	tet_triplets.push_back(SparseMatrixTriplet( 9, 5, v2));	tet_triplets.push_back(SparseMatrixTriplet( 9,  9, v1));
	tet_triplets.push_back(SparseMatrixTriplet(10, 2, v2));	tet_triplets.push_back(SparseMatrixTriplet(10, 6, v2));	tet_triplets.push_back(SparseMatrixTriplet(10, 10, v1));
	tet_triplets.push_back(SparseMatrixTriplet(11, 3, v2));	tet_triplets.push_back(SparseMatrixTriplet(11, 7, v2));	tet_triplets.push_back(SparseMatrixTriplet(11, 11, v1));
	
	m_A_tet.resize(12,12);
	m_A_tet.setFromTriplets( tet_triplets.begin(), tet_triplets.end() );


	////////////////////////////////////////////////////
	// setup B matrices for all constraint types
	////////////////////////////////////////////////////

	m_B_attachment.resize(3,3);
	m_B_attachment = m_A_attachment;

	m_B_spring.resize(6,6);
	m_B_spring = m_A_spring;

	m_B_tet.resize(12,12);
	m_B_tet = m_A_tet;
}

Simulation::~Simulation()
{
	clearConstraints();
}

void Simulation::Reset()
{	
	m_inertia_y.resize(m_mesh->m_system_dimension);
	m_external_force.resize(m_mesh->m_system_dimension);
	
	m_mesh->m_expanded_system_dimension = 0;
	m_mesh->m_expanded_system_dimension_1d = 0;

	setupConstraints();
	CreateLHSMatrix();
	CreateRHSMatrix();
	p_attach.resize(3);
	p_spring.resize(6);
	p_tet.resize(12);
	SetReprefactorFlag();

	for (int i = 0; i < m_mesh->m_vertices_number; i++)
		isColliding.push_back(false);

	m_selected_attachment_constraint = NULL;
	m_step_mode = false;
}

void Simulation::UpdateAnimation(const int fn)
{
	if (m_animation_enable_swinging)
	{
		int swing_num = 0;
		ScalarType swing_step = m_animation_swing_amp / m_animation_swing_half_period;
		EigenVector3 swing_dir(m_animation_swing_dir[0], m_animation_swing_dir[1], m_animation_swing_dir[2]);
		int positive_direction = ((fn/m_animation_swing_half_period)%2)?-1:1;
		for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
		{
			AttachmentConstraint* ac;
			if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
			{
				EigenVector3 new_fixed_point = ac->GetFixedPoint() + swing_dir*swing_step*positive_direction;
				ac->SetFixedPoint(new_fixed_point);
				if (++swing_num >= m_animation_swing_num)
				{
					break;
				}
			}
		}
	}
}

/*********** START - 563 final project primary contribution ***********/


////////////////////////////////////////////////////
// CreateSMatrix()
////////////////////////////////////////////////////
SparseMatrix Simulation::CreateSMatrix(Constraint* c)
{
	SparseMatrix S;
	std::vector<SparseMatrixTriplet> s_triplets;
	AttachmentConstraint* ac;
	if (ac = dynamic_cast<AttachmentConstraint*>(c)) // is attachment constraint
	{
		s_triplets.clear();
		unsigned int m_p0 = ac->GetConstrainedVertexIndex();
		s_triplets.push_back(SparseMatrixTriplet(0, m_p0*3 + 0, 1));
		s_triplets.push_back(SparseMatrixTriplet(1, m_p0*3 + 1, 1));
		s_triplets.push_back(SparseMatrixTriplet(2, m_p0*3 + 2, 1));
		S.resize(3, m_mesh->m_system_dimension);
	}

	SpringConstraint *sc;
	if (sc = dynamic_cast<SpringConstraint*>(c)) // is spring constraint
	{
		s_triplets.clear();
		unsigned int m_p1 = sc->GetConstrainedVertexIndex1();
		s_triplets.push_back(SparseMatrixTriplet(0, m_p1*3 + 0, 1));
		s_triplets.push_back(SparseMatrixTriplet(1, m_p1*3 + 1, 1));
		s_triplets.push_back(SparseMatrixTriplet(2, m_p1*3 + 2, 1));
		unsigned int m_p2 = sc->GetConstrainedVertexIndex2();
		s_triplets.push_back(SparseMatrixTriplet(3, m_p2*3 + 0, 1));
		s_triplets.push_back(SparseMatrixTriplet(4, m_p2*3 + 1, 1));
		s_triplets.push_back(SparseMatrixTriplet(5, m_p2*3 + 2, 1));
		S.resize(6, m_mesh->m_system_dimension);
	}
	
	TetConstraint *tc;
	if (tc = dynamic_cast<TetConstraint*>(c)) // is tetrahedral constraint
	{
		//TODO: all of this
		s_triplets.clear();
		unsigned int m_p1 = tc->GetConstrainedVertexIndex1();
		s_triplets.push_back(SparseMatrixTriplet(0, m_p1*3 + 0, 1));
		s_triplets.push_back(SparseMatrixTriplet(1, m_p1*3 + 1, 1));
		s_triplets.push_back(SparseMatrixTriplet(2, m_p1*3 + 2, 1));
		unsigned int m_p2 = tc->GetConstrainedVertexIndex2();
		s_triplets.push_back(SparseMatrixTriplet(3, m_p2*3 + 0, 1));
		s_triplets.push_back(SparseMatrixTriplet(4, m_p2*3 + 1, 1));
		s_triplets.push_back(SparseMatrixTriplet(5, m_p2*3 + 2, 1));
		unsigned int m_p3 = tc->GetConstrainedVertexIndex3();
		s_triplets.push_back(SparseMatrixTriplet(6, m_p3*3 + 0, 1));
		s_triplets.push_back(SparseMatrixTriplet(7, m_p3*3 + 1, 1));
		s_triplets.push_back(SparseMatrixTriplet(8, m_p3*3 + 2, 1));
		unsigned int m_p4 = tc->GetConstrainedVertexIndex4();
		s_triplets.push_back(SparseMatrixTriplet(9, m_p4*3 + 0, 1));
		s_triplets.push_back(SparseMatrixTriplet(10, m_p4*3 + 1, 1));
		s_triplets.push_back(SparseMatrixTriplet(11, m_p4*3 + 2, 1));
		S.resize(12,m_mesh->m_vertices_number*3);
	}

	S.setFromTriplets(s_triplets.begin(), s_triplets.end());

	return S;
}


////////////////////////////////////////////////////
// CreateLHSMatrix()
////////////////////////////////////////////////////
void Simulation::CreateLHSMatrix()
{
	SparseMatrix M;
	M = m_mesh->m_mass_matrix / (m_h * m_h);

	for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
	{
		SparseMatrix S_i;
		SparseMatrix A_i;
		ScalarType w_i;

		// is attachment constraint
		AttachmentConstraint *ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) {
			w_i = ac->Stiffness();
			A_i = m_A_attachment;
		}

		// is spring constraint
		SpringConstraint *sc;
		if (sc = dynamic_cast<SpringConstraint*>(*c)) {
			w_i = sc->Stiffness();
			A_i = m_A_spring;
		}
	
		// is tetrahedral constraint
		TetConstraint *tc;
		if (tc = dynamic_cast<TetConstraint*>(*c)) {
			w_i = tc->Stiffness();
			A_i = m_A_tet;
		}

		S_i = CreateSMatrix(*c);
		SparseMatrix S_i_transpose = S_i.transpose();
		SparseMatrix A_i_transpose = A_i;
	
		S_i_transpose.applyThisOnTheLeft(A_i_transpose);
		A_i_transpose.applyThisOnTheLeft(A_i);
		A_i.applyThisOnTheLeft(S_i);
			
		M += (w_i * S_i);
	}
	m_llt.compute(M);
}


////////////////////////////////////////////////////
// CreateRHSMatrix()
////////////////////////////////////////////////////
void Simulation::CreateRHSMatrix()
{
	for (int i = 0; i < m_constraints.size(); i++)
	{
		Constraint* c = m_constraints[i];
		SparseMatrix A_i;
		SparseMatrix B_i;
		
		// is attachment constraint
		AttachmentConstraint *ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(c)) {
			A_i = m_A_attachment;
			B_i = m_B_attachment;
		}

		// is spring constraint
		SpringConstraint *sc;
		if (sc = dynamic_cast<SpringConstraint*>(c)) {
			A_i = m_A_spring;
			B_i = m_B_spring;
		}
	
		// is tet constraint
		TetConstraint *tc;
		if (tc = dynamic_cast<TetConstraint*>(c)) {
			A_i = m_A_tet;
			B_i = m_B_tet;
		}

		ScalarType w_i = c->Stiffness();
		SparseMatrix S_i = CreateSMatrix(c);
		SparseMatrix S_i_transpose = S_i.transpose();			
		S_i_transpose.applyThisOnTheLeft(A_i);
		A_i.applyThisOnTheLeft(B_i);
		c->m_RHS = w_i * B_i;
	}
}


////////////////////////////////////////////////////
// Update()
////////////////////////////////////////////////////

void Simulation::Update()
{
	// update inertia term
	calculateInertiaY();

	// update external force
	calculateExternalForce();
	
	// update cloth
	switch ( m_integration_method ) {
		case INTEGRATION_EXPLICIT_EULER:				// TODO
		case INTEGRATION_EXPLICIT_SYMPLECTIC:			// TODO
		case INTEGRATION_IMPLICIT_EULER_BARAFF_WITKIN:	// TODO
		case INTEGRATION_GRADIENT_DESCENT:				// TODO
		case INTEGRATION_NEWTON_DESCENT:				// TODO
		case INTEGRATION_NEWTON_DESCENT_PCG:			// TODO
		case INTEGRATION_PBD:							// TODO
		case INTEGRATION_LOCAL_GLOBAL:					// DONE
		{
			omp_set_num_threads(m_constraints.size());

			AttachmentConstraint* ac;
			SpringConstraint *sc;
			TetConstraint *tc;
					
			VectorX* p_j;
			Constraint* c_j;
			unsigned int tn;
			
			VectorX q_n = m_mesh->m_current_positions;
			VectorX s_n = m_inertia_y + (m_h*m_h)*(m_mesh->m_inv_mass_matrix)*m_external_force;
			VectorX q_n1 = s_n;

			SparseMatrix coeff = m_mesh->m_mass_matrix / (m_h * m_h);
			coeff.applyThisOnTheLeft(s_n);
			
			// LOCAL SOLVE STEP
			for (int i = 0; i < m_iterations_per_frame; i++)
			{
				simpleTimer localTimer;

				VectorX b = s_n;

				EigenVector3 current_vector;
				ScalarType current_stretch;
				int constraintType;

				int num_parallel_loops = ceil( m_constraints.size() / (float) omp_get_max_threads() );
				#pragma omp for
				for (int j = 0; j < num_parallel_loops; j++)
				{
					#pragma omp parallel default(shared) private(c_j, p_j, tn, sc, ac, tc, current_stretch, current_vector)
					{
						tn = omp_get_thread_num() + j*omp_get_max_threads();
						if (tn < m_constraints.size())
						{
							c_j = m_constraints[tn];
							constraintType = c_j->constraintType;

							#pragma omp critical
							{
								if (constraintType == SPRING) // is spring constraint
								{
									sc = (SpringConstraint *) c_j;
									current_vector = q_n1.block_vector(sc->GetConstrainedVertexIndex1()) - q_n1.block_vector(sc->GetConstrainedVertexIndex2());
									current_stretch = current_vector.norm() - sc->GetRestLength();
									current_vector = (current_stretch / 2.0) * current_vector.normalized();
									
									p_j = &p_spring;
									p_j->block_vector(0) = q_n1.block_vector(sc->GetConstrainedVertexIndex1()) - current_vector;
									p_j->block_vector(1) = q_n1.block_vector(sc->GetConstrainedVertexIndex2()) + current_vector;
								}

								else if (constraintType == ATTACHMENT) // is attachment constraint
								{
									ac = (AttachmentConstraint *) c_j;
									p_j = &p_attach;
									p_j->block_vector(0) = ac->GetFixedPoint();
								}

								else if (constraintType == TET) // is tetrahedral constraint
								{
									tc = ( TetConstraint* )c_j;

									// build VectorX with deformed tetrahedron's current vertex positions
									VectorX tet_verts_current;
									tet_verts_current.resize( 12 );
									tet_verts_current.block_vector( 0 ) = q_n1.block_vector( tc->GetConstrainedVertexIndex1() );
									tet_verts_current.block_vector( 1 ) = q_n1.block_vector( tc->GetConstrainedVertexIndex2() );
									tet_verts_current.block_vector( 2 ) = q_n1.block_vector( tc->GetConstrainedVertexIndex3() );
									tet_verts_current.block_vector( 3 ) = q_n1.block_vector( tc->GetConstrainedVertexIndex4() );

									// call computeVolumePreservingVertexPositions() to compute deformed tetrahedron's new vertex positions that preserve volume
									VectorX tet_verts_new;
									tet_verts_new.resize( 12 );
									tc->computeVolumePreservingVertexPositions( tet_verts_new, tet_verts_current );

									// set deformed tetrahedron's new vertex positions that preserve volume
									p_j = &p_tet;
									p_j->block_vector( 0 ) = tet_verts_new.block_vector( 0 );
									p_j->block_vector( 1 ) = tet_verts_new.block_vector( 1 );
									p_j->block_vector( 2 ) = tet_verts_new.block_vector( 2 );
									p_j->block_vector( 3 ) = tet_verts_new.block_vector( 3 );
								}

								c_j->m_RHS.applyThisOnTheLeft(*p_j);
								b += *p_j;
							}
						}
					}
				}			
				// GLOBAL SOLVE STEP
				q_n1 = m_llt.solve(b);
			}

			VectorX v_n1 = (q_n1 - q_n)/m_h;
			m_mesh->m_current_positions = q_n1;
			m_mesh->m_current_velocities = v_n1;
		
			break;
		}
	}
	
	// Add collision detection here using pos_next;
	VectorX penetration = collisionDetection(m_mesh->m_current_positions);
	m_mesh->m_current_positions -= penetration;

	// update velocity and damp
	dampVelocity();
}

void Simulation::DrawConstraints(const VBO& vbos)
{
	for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
	{
		(*it)->Draw(vbos);
	}
}

ScalarType Simulation::TryToSelectAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
	ScalarType ray_point_dist;
	ScalarType min_dist = 100.0;
	AttachmentConstraint* best_candidate = NULL;

	bool current_state_on = false;
	for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
	{
		AttachmentConstraint* ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
		{
			ray_point_dist = ((ac->GetFixedPoint()-p0).cross(dir)).norm();
			if (ray_point_dist < min_dist)
			{
				min_dist = ray_point_dist;
				best_candidate = ac;
			}
		}
	}
	// exit if no one fits
	if (min_dist > DEFAULT_SELECTION_RADIUS)
	{
		UnselectAttachmentConstraint();

		return -1;
	}
	else
	{
		SelectAtttachmentConstraint(best_candidate);
		EigenVector3 fixed_point_temp = m_mesh->m_current_positions.block_vector(m_selected_attachment_constraint->GetConstrainedVertexIndex());

		return (fixed_point_temp-p0).dot(dir); // this is m_cached_projection_plane_distance
	}
}

bool Simulation::TryToToggleAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
	EigenVector3 p1;

	ScalarType ray_point_dist;
	ScalarType min_dist = 100.0;
	unsigned int best_candidate = 0;
	// first pass: choose nearest point
	for (unsigned int i = 0; i != m_mesh->m_vertices_number; i++)
	{
		p1 = m_mesh->m_current_positions.block_vector(i);

		ray_point_dist = ((p1-p0).cross(dir)).norm();
		if (ray_point_dist < min_dist)
		{
			min_dist = ray_point_dist;
			best_candidate = i;
		}
	}
	for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
	{
		AttachmentConstraint* ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
		{
			ray_point_dist = ((ac->GetFixedPoint()-p0).cross(dir)).norm();
			if (ray_point_dist < min_dist)
			{
				min_dist = ray_point_dist;
				best_candidate = ac->GetConstrainedVertexIndex();
			}
		}
	}
	// exit if no one fits
	if (min_dist > DEFAULT_SELECTION_RADIUS)
	{
		return false;
	}
	// second pass: toggle that point's fixed position constraint
	bool current_state_on = false;
	for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
	{
		AttachmentConstraint* ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
		{
			if (ac->GetConstrainedVertexIndex() == best_candidate)
			{
				current_state_on = true;
				m_constraints.erase(c);
				delete ac;
				m_mesh->m_expanded_system_dimension-=3;
				m_mesh->m_expanded_system_dimension_1d-=1;
				CreateLHSMatrix();
				CreateRHSMatrix();
				break;
			}
		}
	}
	if (!current_state_on)
	{
		AddAttachmentConstraint(best_candidate);
	}

	return true;
}

void Simulation::SelectAtttachmentConstraint(AttachmentConstraint* ac)
{
	m_selected_attachment_constraint = ac;
	m_selected_attachment_constraint->Select();
}

void Simulation::UnselectAttachmentConstraint()
{
	if (m_selected_attachment_constraint)
	{
		m_selected_attachment_constraint->UnSelect();
	}
	m_selected_attachment_constraint = NULL;
}

void Simulation::AddAttachmentConstraint(unsigned int vertex_index)
{
	AttachmentConstraint* ac = new AttachmentConstraint(&m_stiffness_attachment, vertex_index, m_mesh->m_current_positions.block_vector(vertex_index));
	m_constraints.push_back(ac);
	m_mesh->m_expanded_system_dimension+=3;
	m_mesh->m_expanded_system_dimension_1d+=1;
	CreateLHSMatrix();
	CreateRHSMatrix();
}

void Simulation::MoveSelectedAttachmentConstraintTo(const EigenVector3& target)
{
	if (m_selected_attachment_constraint)
		m_selected_attachment_constraint->SetFixedPoint(target);
}

void Simulation::SaveAttachmentConstraint(const char* filename)
{
	std::ofstream outfile;
	outfile.open(filename, std::ifstream::out);
	if (outfile.is_open())
	{
		int existing_vertices = 0;
		for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
		{
			(*c)->WriteToFileOBJHead(outfile);
		}
		outfile << std::endl;
		for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
		{
			(*c)->WriteToFileOBJ(outfile, existing_vertices);
		}

		outfile.close();
	}
}

void Simulation::LoadAttachmentConstraint(const char* filename)
{
	// clear current attachement constraints
	for (std::vector<Constraint*>::iterator& c = m_constraints.begin(); c != m_constraints.end(); )
	{
		AttachmentConstraint* ac;
		if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
		{
			c = m_constraints.erase(c);
			delete ac;
			m_mesh->m_expanded_system_dimension-=3;
			m_mesh->m_expanded_system_dimension_1d-=1;
		}
		else
		{
			c++;
		}
	}

	// read from file
	std::ifstream infile;
	infile.open(filename, std::ifstream::in);
	char ignore[256];
	if (infile.is_open())
	{
		while(!infile.eof())
		{
			int id;
			EigenVector3 p;
			if (infile >> ignore >> id >> p[0] >> p[1] >> p[2])
			{
				if (strcmp(ignore, "v") == 0)
					break;
				AttachmentConstraint* ac = new AttachmentConstraint(&m_stiffness_attachment, id, p);
				m_constraints.push_back(ac);
				m_mesh->m_expanded_system_dimension+=3;
				m_mesh->m_expanded_system_dimension_1d+=1;
			}
			else
				break;
		}

		infile.close();
	}

}

void Simulation::clearConstraints()
{
	for (unsigned int i = 0; i < m_constraints.size(); ++i)
	{
		delete m_constraints[i];
	}
	m_constraints.clear();
}

void Simulation::setupConstraints()
{
	clearConstraints();

	switch(m_mesh->m_mesh_type)
	{
	case MESH_TYPE_ROPE:
	{
		EigenVector3 pR1, pR2;
		for (std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
		{
			pR1 = m_mesh->m_current_positions.block_vector(e->m_v1);
			pR2 = m_mesh->m_current_positions.block_vector(e->m_v2);
			SpringConstraint* c = new SpringConstraint(&m_stiffness_stretch, e->m_v1, e->m_v2, (pR1 - pR2).norm());
			m_constraints.push_back(c);
			m_mesh->m_expanded_system_dimension += 6;
			m_mesh->m_expanded_system_dimension_1d += 2;
		}

		// Setting up the bending constraints
		for (std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
		{
			pR1 = m_mesh->m_current_positions.block_vector(e->m_v1);
			pR2 = m_mesh->m_current_positions.block_vector(e->m_v2);
			SpringConstraint* c = new SpringConstraint(&m_stiffness_bending, e->m_v1, e->m_v2, (pR1 - pR2).norm());
			m_constraints.push_back(c);
			m_mesh->m_expanded_system_dimension += 6;
			m_mesh->m_expanded_system_dimension_1d += 2;
		}

		AddAttachmentConstraint(0);
		// Add edge to mesh.
	}
		break;
	case MESH_TYPE_CLOTH:
		// procedurally generate constraints including to attachment constraints
		{
			// generate stretch constraints. assign a stretch constraint for each edge.
			EigenVector3 p1, p2;
			for(std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
			{
				p1 = m_mesh->m_current_positions.block_vector(e->m_v1);
				p2 = m_mesh->m_current_positions.block_vector(e->m_v2);
				SpringConstraint* c = new SpringConstraint(&m_stiffness_stretch, e->m_v1, e->m_v2, (p1-p2).norm());
				m_constraints.push_back(c);
				m_mesh->m_expanded_system_dimension+=6;
				m_mesh->m_expanded_system_dimension_1d+=2;
			}

			// generate bending constraints. naive way
			unsigned int i, k;
			for(i = 0; i < m_mesh->m_dim[0]; ++i)
			{
				for(k = 0; k < m_mesh->m_dim[1]; ++k)
				{
					unsigned int index_self = m_mesh->m_dim[1] * i + k;
					p1 = m_mesh->m_current_positions.block_vector(index_self);
					if (i+2 < m_mesh->m_dim[0])
					{
						unsigned int index_row_1 = m_mesh->m_dim[1] * (i + 2) + k;
						p2 = m_mesh->m_current_positions.block_vector(index_row_1);
						SpringConstraint* c = new SpringConstraint(&m_stiffness_bending, index_self, index_row_1, (p1-p2).norm());
						m_constraints.push_back(c);
						m_mesh->m_expanded_system_dimension+=6;
						m_mesh->m_expanded_system_dimension_1d+=2;
					}
					if (k+2 < m_mesh->m_dim[1])
					{
						unsigned int index_column_1 = m_mesh->m_dim[1] * i + k + 2;
						p2 = m_mesh->m_current_positions.block_vector(index_column_1);
						SpringConstraint* c = new SpringConstraint(&m_stiffness_bending, index_self, index_column_1, (p1-p2).norm());
						m_constraints.push_back(c);
						m_mesh->m_expanded_system_dimension+=6;
						m_mesh->m_expanded_system_dimension_1d+=2;
					}
				}
			}

			// generating attachment constraints.
			AddAttachmentConstraint(0);
			 
		}
		break;
	case MESH_TYPE_TET:
		{
			VectorX& x = m_mesh->m_current_positions;
			TetMesh* tet_mesh = dynamic_cast<TetMesh*>(m_mesh);
			for (unsigned int i = 0; i < tet_mesh->m_loaded_mesh->m_tets.size(); ++i)
			{
				MeshLoader::Tet& tet = tet_mesh->m_loaded_mesh->m_tets[i];
				TetConstraint *c = new TetConstraint(&m_stiffness_stretch, tet.id1, tet.id2, tet.id3, tet.id4, x);
				m_constraints.push_back(c);
				m_mesh->m_expanded_system_dimension+=9;
				m_mesh->m_expanded_system_dimension_1d+=3;
			}
		}
		break;
	}
}

void Simulation::dampVelocity()
{
	// TODO: post processing damping
	m_mesh->m_current_velocities *= (1.0 - m_damping_coefficient);
}

void Simulation::calculateInertiaY()
{
	m_inertia_y = m_mesh->m_current_positions + m_mesh->m_current_velocities * m_h;
}

void Simulation::calculateExternalForce()
{
	m_external_force.resize(m_mesh->m_system_dimension);
	m_external_force.setZero();

	// gravity
	for (unsigned int i = 0; i < m_mesh->m_vertices_number; ++i)
	{
		m_external_force[3*i+1] += -m_gravity_constant;
	}

	m_external_force = m_mesh->m_mass_matrix * m_external_force;
}

VectorX Simulation::collisionDetection(const VectorX x)
{
	// Naive implementation of collision detection
	VectorX penetration(m_mesh->m_system_dimension);
	penetration.setZero();
	EigenVector3 normal;
	ScalarType dist;

	EigenVector3 vn, vt, vel;

	float friction = 0.98f;
	float restitution = 0.4f;


	for (unsigned int i = 0; i < m_mesh->m_vertices_number; ++i)
	{
		EigenVector3 xi = x.block_vector(i);

		if (m_scene->StaticIntersectionTest(xi, normal, dist))
		{
			penetration.block_vector( i ) +=  dist * normal;

			isColliding[i] = true;

			// TODO: update velocity

			vel = m_mesh->m_current_velocities.block_vector( i );

			// compute component of velocity parallel to normal at collision surface point
			// glm::dot( vel, normal ) = signed length of projection of vel onto any line parallel to normal
			// divide by glm::dot( normal, normal ) to remove any scaling caused by normal's magnitude
			// multiply by normal to convert previously computed scalar into a vector with direction parallel to normal
			vn = vel.dot( normal ) / normal.dot( normal ) * normal;

			// compute component of velocity perpendicular to normal at collision surface point
			vt = vel - vn;

			// damp velocity components
			vn = vn * restitution;
			vt = vt * friction;

			// reflect velocity component parallel to collision normal
			vn = vn * -1.0f;

			// set particle velocity to dampened velocity
			m_mesh->m_current_velocities.block_vector( i ) = vn + vt;
		}

		else isColliding[i] = false;
	}

	return penetration;
}
