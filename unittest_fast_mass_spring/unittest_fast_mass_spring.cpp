// unittest_fast_mass_spring.cpp : Defines the entry point for the console application.
//
#include "gtest/gtest.h"
#include "mesh.h"
#include "converge.cu"

int _tmain(int argc)
{
	return 0;
}

TEST(EnsureRopeInitialisesWith100Vertices, testGetNumberOfVertices)
{
	RopeMesh* mesh = new RopeMesh(100);
	ASSERT_EQ(100, mesh->GetNumberOfVertices());
}



TEST(EnsureRopeInitalisesWithTrue, testRopeInit)
{
	RopeMesh* mesh = new RopeMesh(100);
	ASSERT_TRUE(mesh->Init());
}


TEST(EnsureDouble3ToDoubleWorks, double3ToDouble) {
	double3 a = make_double3(10, 10, 10);
	double real[3] = { 10, 10, 10 };
	ASSERT_EQ(real, double3ToDouble(a));
}
