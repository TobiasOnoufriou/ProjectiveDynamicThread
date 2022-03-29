ProjectionBasedDynamics
=======================
Repository is forked from: https://github.com/rarietta/Projective-Dynamics

Implementation for final year project:
* Includes partial implementation of CUDA local solver
* Tetrahedral meshes generated using Tetgen in filepath: fast_mass_spring/tet_rope_mesh 
* 1D/2D line in 3D space to act as the thread since volumetric mesh would cause performance issues.


Files contributed to: 
CudaConstraint.cpp, CudaConstraint.h, converge.cu, converge.cuh, simpleTimer.cpp (Change return type for averageComputation), simpleTimer.h and simulation.h, simulation.cpp and anttweakbar_wrapper.cpp (for UI to change the amount of joints of the thread)


Source files are in fast_mass_spring/source

test files are in unittest_fast_mass_spring/unittest_fast_mass_spring.cpp
