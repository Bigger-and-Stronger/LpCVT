//
// This file is derived from LpCVT/main.cpp in the original code
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-1.
//

//===================================================================================================
//
// Algebra tests
//
//  Computes :
//     F_{L_p} in the surface case
//     F_{L_p} in the volume case
//
//===================================================================================================

#ifndef GEEX_TEST_ALGEBRA_H
#define GEEX_TEST_ALGEBRA_H

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/clipped_VD.h>
#include <LpCVT/algebra/F_Lp.h>
#include <LpCVT/common/line_stream.h>
#include <fstream>

namespace Geex {
	/** Used by get_combinatorics() in volume mode **/
	class MemorizeIndices {
	public:
		MemorizeIndices(std::vector<int>& I_in, std::vector<vec3>& C_in) : I(I_in), C(C_in) {
			std::vector<int>().swap(I);
			std::vector<vec3>().swap(C);
		}

		void operator() (const unsigned int i, int j,
			const VertexEdge& v1, const VertexEdge& v2, const VertexEdge& v3) const;
	private:
		std::vector<int>& I;
		std::vector<vec3>& C;
	};

	/** Used by get_combinatorics() in surface mode **/
	class MemorizeIndicesAndFacets {
	public:
		MemorizeIndicesAndFacets(const RestrictedVoronoiDiagram& RVD_in,
			std::vector<int>& I_in, std::vector<vec3>& C_in, std::vector<int>& F_in) : RVD(RVD_in), I(I_in), C(C_in), F(F_in) {
			std::vector<int>().swap(I);
			std::vector<vec3>().swap(C);
			std::vector<int>().swap(F);
		}

		void operator() (const unsigned int i,
			const VertexEdge& v1, const VertexEdge& v2, const VertexEdge& v3) const;
	private:
		const RestrictedVoronoiDiagram& RVD;
		std::vector<int>& I;
		std::vector<vec3>& C;
		std::vector<int>& F;
	};

	/** Gets the combinatorics of the integration simplices,
	 * i.e. 10 integers per integration simplex.
	 * (see Section 3.1 in the paper)
	 * Returns also the array of C vertices (three per integration simplex).
	 * Since they are easy to get during the combinatorial phase, they are
	 * computed here and kept for the algebraic phase.
	 *
	 * In 2D mode (volume = false), returns also the array F.
	 *   F[i] indicates the facet that contains the i-th integration simplex. **/
	void get_combinatorics(
		Mesh* M, const std::vector<vec3>& pts,
		std::vector<int>& I, std::vector<vec3>& C, std::vector<int>& F, bool volume);

	/** Computes F_{L_p} and its gradient. **/
	void compute_F_g(
		Mesh* m, const std::vector<vec3>& pts,
		unsigned int p, bool volume);

	void test_algebra(const std::string& mesh_filename, const std::string& pts_filename);
}

#endif //GEEX_PROCESSOR_H
