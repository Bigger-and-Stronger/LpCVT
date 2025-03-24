//
// This file is derived from LpCVT/main.cpp in the original code from https://xueyuhanlang.github.io
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-1.
//

#include "LpCVT/others/Geex_test_algebra.h"
#include "LpCVT/others/Geex_test_combinatorics.h"
#include "LpCVT/others/macro.h"

namespace Geex {
	/** Used by get_combinatorics() in volume mode **/
	class MemorizeIndices {
	public:
		MemorizeIndices(std::vector<int>& I_in, std::vector<vec3>& C_in) : I(I_in), C(C_in) {
			std::vector<int>().swap(I);
			std::vector<vec3>().swap(C);
		}

		void operator() (
			const unsigned int i, int j,
			const VertexEdge& v1, const VertexEdge& v2, const VertexEdge& v3
			) const {
			I.push_back(i);
			I.push_back(v1.sym[2]); I.push_back(v1.sym[1]); I.push_back(v1.sym[0]);
			I.push_back(v2.sym[2]); I.push_back(v2.sym[1]); I.push_back(v2.sym[0]);
			I.push_back(v3.sym[2]); I.push_back(v3.sym[1]); I.push_back(v3.sym[0]);
			C.push_back(v1); C.push_back(v2); C.push_back(v3);
		}
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

		void operator() (
			const unsigned int i,
			const VertexEdge& v1, const VertexEdge& v2, const VertexEdge& v3
			) const {
			I.push_back(i);
			I.push_back(v1.sym[2]); I.push_back(v1.sym[1]); I.push_back(v1.sym[0]);
			I.push_back(v2.sym[2]); I.push_back(v2.sym[1]); I.push_back(v2.sym[0]);
			I.push_back(v3.sym[2]); I.push_back(v3.sym[1]); I.push_back(v3.sym[0]);
			F.push_back(RVD.current_facet());
			C.push_back(v1); C.push_back(v2); C.push_back(v3);
		}
	private:
		const RestrictedVoronoiDiagram& RVD;
		std::vector<int>& I;
		std::vector<vec3>& C;
		std::vector<int>& F;
	};

	void get_combinatorics(
		Mesh* M, const std::vector<vec3>& pts,
		std::vector<int>& I, std::vector<vec3>& C, std::vector<int>& F,
		const bool volume
	) {
		Delaunay* delaunay = Delaunay::create("CGAL");
		delaunay->set_vertices(pts);

		if (volume) {
			ClippedVoronoiDiagram CVD(delaunay, M);
			CVD.for_each_triangle(MemorizeIndices(I, C));
		}
		else {
			RestrictedVoronoiDiagram RVD(delaunay, M);
			RVD.set_symbolic(true);
			RVD.for_each_triangle(MemorizeIndicesAndFacets(RVD, I, C, F));
		}
		delete delaunay;
	}

	void compute_F_g(
		Mesh* m, const std::vector<vec3>& pts,
		const unsigned int p,
		const bool volume
		) {
		VERBOSE("nb pts = " << pts.size() << " nb facets = " << m->nb_facets());
		std::vector<int> I;
		std::vector<vec3> C;
		std::vector<int> F;
		get_combinatorics(m, pts, I, C, F, volume);
		const unsigned int nb_integration_simplices = static_cast<unsigned int>(I.size()) / 10;
		std::vector<mat3> M(nb_integration_simplices);
		for (auto & i : M) {
			i.load_identity();
			// or replace with anisotropy field
			//   In 2D: use F[i] to retrieve the index of the facet that contains
			//      the current integration simplex (and access an array of per-facet anisotropy).
			//   In 3D: use geometric search from the centroid of the current
			//      integration simplex.
		}
		std::vector<plane3> Q(m->nb_facets());
		for (unsigned int i = 0; i < m->nb_facets(); ++i)
			Q[i] = m->facet_plane(i);

		std::vector<double> g(pts.size() * 3);
		const double f = compute_F_Lp(volume, p, m, I, C, pts, Q, M, g);
		double gnorm = 0.0;
		for (const double i : g)
			gnorm += i * i;

		gnorm = ::sqrt(gnorm);
		std::cout.precision(16);

		VERBOSE(
			(volume ? "volume " : "surface ")
			<< "F_L" << p << ":"
			<< "f=" << std::scientific << f << "  g=" << gnorm
			);
	}

    void test_algebra(
    	const std::string& mesh_filename, const std::string& pts_filename,
    	const int p
    	) {
        Mesh M;
        const unsigned int nb_borders = M.load(mesh_filename);
        std::vector<vec3> pts;
        if (!load_pts(pts_filename, pts))
            return;

        if (nb_borders == 0) {
            VERBOSE("========== volume LpCVT test ==========");
            compute_F_g(&M, pts, p, true);
        }

        VERBOSE("========== surface LpCVT test ==========");
        compute_F_g(&M, pts, p, false);
    }
}