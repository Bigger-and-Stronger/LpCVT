//
// This file is derived from LpCVT/main.cpp in the original code
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-1.
//

#include "LpCVT/others/Geex_test_algebra.h"
#include "LpCVT/others/Geex_test_combinatorics.h"
#include "LpCVT/others/macro.h"

namespace Geex {
	void MemorizeIndices::operator() (const unsigned int i, int j,
			const VertexEdge& v1, const VertexEdge& v2, const VertexEdge& v3
			) const {
		I.push_back(i);
		I.push_back(v1.sym[2]); I.push_back(v1.sym[1]); I.push_back(v1.sym[0]);
		I.push_back(v2.sym[2]); I.push_back(v2.sym[1]); I.push_back(v2.sym[0]);
		I.push_back(v3.sym[2]); I.push_back(v3.sym[1]); I.push_back(v3.sym[0]);
		C.push_back(v1); C.push_back(v2); C.push_back(v3);
	}

	void MemorizeIndicesAndFacets::operator() (const unsigned int i,
			const VertexEdge& v1, const VertexEdge& v2, const VertexEdge& v3
			) const {
		I.push_back(i);
		I.push_back(v1.sym[2]); I.push_back(v1.sym[1]); I.push_back(v1.sym[0]);
		I.push_back(v2.sym[2]); I.push_back(v2.sym[1]); I.push_back(v2.sym[0]);
		I.push_back(v3.sym[2]); I.push_back(v3.sym[1]); I.push_back(v3.sym[0]);
		F.push_back(RVD.current_facet());
		C.push_back(v1); C.push_back(v2); C.push_back(v3);
	}

	void get_combinatorics(
		Mesh* M, const std::vector<vec3>& pts,
		std::vector<int>& I, std::vector<vec3>& C, std::vector<int>& F, const bool volume
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
		unsigned int p, const bool volume
		) {
		VERBOSE("nb pts = " << pts.size() << "   nb facets = " << m->nb_facets());
		std::vector<int> I;
		std::vector<vec3> C;
		std::vector<int> F;
		get_combinatorics(m, pts, I, C, F, volume);
		const unsigned int nb_integration_simplices = static_cast<unsigned int>(I.size()) / 10;
		std::vector<mat3> M(nb_integration_simplices);
		for (unsigned int i = 0; i < M.size(); ++i) {
			M[i].load_identity();
			// or replace with anisotropy field
			//   In 2D: use F[i] to retreive the index of the facet that contains
			//      the current integration simplex (and access an array of per-facet anisotropy).
			//   In 3D: use geometric search from the centroid of the current
			//      integration simplex.
		}
		std::vector<plane3> Q(m->nb_facets());
		for (unsigned int i = 0; i < m->nb_facets(); ++i)
			Q[i] = m->facet_plane(i);

		std::vector<double> g(pts.size() * 3);
		double f = compute_F_Lp(volume, p, m, I, C, pts, Q, M, g);
		double gnorm = 0.0;
		for (unsigned int i = 0; i < g.size(); ++i)
			gnorm += g[i] * g[i];

		gnorm = ::sqrt(gnorm);
		std::cout.precision(16);
		VERBOSE((volume ? "volume " : "surface ")
			<< "F_L" << p << ":"
			<< "f=" << std::scientific << f << "  g=" << gnorm);
	}

    void test_algebra(
    	const std::string& mesh_filename, const std::string& pts_filename
    	) {
        Mesh M;
        const unsigned int nb_borders = M.load(mesh_filename);
        std::vector<vec3> pts;
        if (!load_pts(pts_filename, pts))
            return;

        if (nb_borders == 0) {
            VERBOSE("========== volume LpCVT test ======");
            compute_F_g(&M, pts, 4, true);
        }
        VERBOSE("========== surface LpCVT test ======");
        compute_F_g(&M, pts, 4, false);
    }
}