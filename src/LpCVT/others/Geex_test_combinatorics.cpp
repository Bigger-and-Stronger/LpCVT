//
// This file is derived from LpCVT/main.cpp in the original code from https://xueyuhanlang.github.io
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-1.
//

#include "LpCVT/others/Geex_test_combinatorics.h"
#include "LpCVT/others/IO.h"
#include "LpCVT/others/macro.h"

#include <LpCVT/common/line_stream.h>

namespace Geex {
	// ------------------------------ RDT ------------------------------
	/** Used by save_RDT(). **/
	class SavePrimalTriangle {
	public:
		explicit SavePrimalTriangle(std::ofstream& out): out_(&out) {}

		void operator()(const unsigned int i, const unsigned int j, const unsigned int k
			) const {
			(*out_) << "f " << i + 1 << " " << j + 1 << " " << k + 1 << std::endl;
		}
	private:
		std::ofstream* out_;
	};

	void save_RDT(
		RestrictedVoronoiDiagram& RVD, const std::string& filename
		) {
		VERBOSE("Computing and saving RDT to " << filename);
		std::ofstream out(filename.c_str());
		for (unsigned int i = 0; i < RVD.delaunay()->nb_vertices(); ++i)
			out << "v " << RVD.delaunay()->vertex(i) << std::endl;

		RVD.for_each_primal_triangle(SavePrimalTriangle(out));
		out.close();
	}

	// ------------------------------ RVD ------------------------------
	/** used by save_RVD() **/
	class SaveRVDFacets {
	public:
		explicit SaveRVDFacets(std::ostream& out) : out_(out), cur_v_(1), cur_f_(1) {
			out << "# attribute chart facet integer" << std::endl;
		}
		void operator()(const unsigned int iv, Mesh* M
			) const {
			for (unsigned int f = 0; f < M->nb_facets(); f++) {
				for (unsigned int i = M->facet_begin(f); i < M->facet_end(f); ++i) {
					const vec3& v = M->vertex(i);
					out_ << "v " << v << std::endl;
				}
				out_ << "f ";
				for (unsigned int i = M->facet_begin(f); i < M->facet_end(f); ++i) {
					out_ << cur_v_++ << " ";
				}
				out_ << std::endl;
				out_ << "# attrs f " << cur_f_++ << " " << iv << std::endl;
			}
		}
	private:
		std::ostream& out_;
		mutable unsigned int cur_v_; // vertex cnt
		mutable unsigned int cur_f_; // face cnt
	};

	void save_RVD(
		RestrictedVoronoiDiagram& RVD, const std::string& filename
		) {
		std::ofstream out(filename.c_str());
		if (!out) {
			WARNING("could not open " << filename);
			return;
		}

		VERBOSE("Computing and saving RVD to " << filename);
		const bool sym_backup = RVD.symbolic();
		RVD.set_symbolic(true);
		RVD.for_each_facet(SaveRVDFacets(out));
		RVD.set_symbolic(sym_backup);
	}

	// ------------------------------ CVD ------------------------------

	void SaveClippedVDFacets::operator()(
		const unsigned int i, int j, const vec3& p1, const vec3& p2, const vec3& p3
		) const {
		const vec3 x0 = delaunay_->vertex(i);
		out_ << "v " << x0 + shrink_ * (p1 - x0) << std::endl;
		out_ << "v " << x0 + shrink_ * (p2 - x0) << std::endl;
		out_ << "v " << x0 + shrink_ * (p3 - x0) << std::endl;
		out_ << "f " << cur_ << " " << cur_ + 1 << " " << cur_ + 2 << std::endl;
		cur_ += 3;
		out_ << "# attrs f " << ((cur_ - 1) / 3) << " " << i << std::endl;
	}

	void save_clippedVD(
		ClippedVoronoiDiagram& CVD, const std::string& filename, double shrink
		) {
		std::ofstream out(filename.c_str());
		if (!out) {
			WARNING("could not open " << filename);
			return;
		}
		VERBOSE("Computing and saving clipped VD to " << filename);
		CVD.for_each_triangle(SaveClippedVDFacets(CVD.delaunay(), out, shrink));
	}



	void test_combinatorics(
		const std::string& mesh_filename, const std::string& pts_filename,
		const std::string& output_path
		) {
		Mesh M;
		const unsigned int nb_borders = M.load(mesh_filename);
		std::vector<vec3> pts;
		if (!load_pts(pts_filename, pts))
			return;

		Delaunay* delaunay = Delaunay::create("CGAL");
		RestrictedVoronoiDiagram RVD(delaunay, &M);
		ClippedVoronoiDiagram CVD(delaunay, &M);

		delaunay->set_vertices(pts);
		save_RVD(RVD, output_path + "rvd.obj");
		save_RDT(RVD, output_path + "rdt.obj");
		if (nb_borders == 0)
			save_clippedVD(CVD, output_path + "cvd.obj", 0.7);

		delete delaunay;
	}
}