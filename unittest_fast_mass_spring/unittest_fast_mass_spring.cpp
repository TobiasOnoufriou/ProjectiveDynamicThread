// unittest_fast_mass_spring.cpp : Defines the entry point for the console application.
// Code below created by Tobias Onoufriou
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

TEST(EnsureAjoinVectorArrayWorks, aJoinVectorArray) {
	double a[3] = { 10, 3 ,4 };
	double b[3] = { 5, 4,2 };
	double c[6] = { 10, 3, 4, 5, 4, 2 };
	double* out;

	ajoinVectorArray(a, b, out);

	EXPECT_TRUE(0 == std::memcmp(c, out, sizeof(c)));
}

TEST(EnsureNormaliseWorks, normalise) {
	double3 a = make_double3(3, 4, 2);
	
	double3 expected = make_double3(0.5570860145311556, 0.7427813527082074, 0.3713906763541037);
	ASSERT_EQ(expected, normalise(a));
}