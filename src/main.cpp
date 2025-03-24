// Quick reference:
//
// Mesh: stores a Piecewise Linear Complex (+ functions to read it from a file in .obj format)
// Delaunay: abstract interface to Delaunay triangulation (+ Delaunay_CGAL: implementation with CGAL)
// RestrictedVoronoiDiagram: Given a Mesh and a Delaunay, computes the intersection between a Voronoi diagram and a Mesh.
// ClippedVoronoiDiagram: computes the intersection between a Voronoi diagram and the interior of a closed Mesh
//
// RestrictedVoronoiDiagram and ClippedVoronoiDiagram are accessed through the for_each_triangle(do_it) function,
// that does the actual computation, and calls the user-provided function do_it(i,j,k,l) for each integration
// simplex. Note that do_it can be a function or an object that overloards operator().
// This mechanism can be used for both computing F-Lp and displaying the clipped Voronoi cells as in figures 4,5,6
// in the paper.

// This code has been modified by Canjia Huang <canjia7@gmail.com> on 25-3-1.

#include "LpCVT/others/Geex_test_combinatorics.h"
#include "LpCVT/others/Geex_test_algebra.h"
#include "LpCVT/others/functional.h"
#include "LpCVT/others/IO.h"
#include "LpCVT/others/optimizer.h"
#include "LpCVT/others/macro.h"
#include "LpCVT/others/LBFGSpp_solver.h"
#include <LpCVT/common/types.h>
#include <LpCVT/combinatorics/mesh.h>
#include "LpCVT/external/CLI11/CLI11.hpp"

#include <filesystem>

int main(int argc, char** argv) {
	std::string filename_mesh;
	std::string filename_pts;
	bool volume = false;
	unsigned int p = 2;
	unsigned int max_iter = 100;
	unsigned int sample_num = 10000;

	CLI::App app{"LpCVT"};
	argv = app.ensure_utf8(argv);

	app.add_option("-p,-P", p, "Set L_p, p must be even, default: 2");
	app.add_flag("-v,-V", volume, "Calculate in volume, default: fasle");
	app.add_option("-m,-M", max_iter, "Maximum number of iterations, default: 100");
	app.add_option("-s,-S", sample_num, "The number of randomly sampled points on the input mesh, default: 1e4");
	app.add_option("meshPath", filename_mesh, "Reference mesh path")->check(CLI::ExistingFile)->required();
	app.add_option("ptsPath", filename_pts, "Pts path")->check(CLI::ExistingFile);

	CLI11_PARSE(app, argc, argv);

	/** Process **/
	const std::filesystem::path absolute_path = std::filesystem::absolute(filename_mesh);
	std::filesystem::path parent_dir = absolute_path.parent_path();
	const std::string output_path = parent_dir.string() + "/";

	if (p % 2 != 0) {
		std::cout << "error: p must be even" << std::endl;
		return 1;
	}

	Geex::Mesh M;
	M.load(filename_mesh);

	if (filename_pts.empty()) {
		std::vector<Geex::vec3> pts;
		Geex::sample_mesh_pts(M, pts, sample_num);
		Geex::write_pts(output_path + "sample.pts", pts);
		filename_pts = output_path + "sample.pts";
	}

	std::vector<Geex::vec3> pts;
	if (!Geex::load_pts(filename_pts, pts))
		return 1;

	Geex::optimize(LBFGSppSolver::LBFGSpp_solver(&M, pts, volume, p, max_iter));

	Geex::write_pts(output_path + "res_sample.pts", pts);
	filename_pts = output_path + "res_sample.pts";

	std::cout << "============= geometry -> combinatorics test ==========" << std::endl;
	Geex::test_combinatorics(filename_mesh, filename_pts, output_path);
	std::cout << std::endl;
	std::cout << "============= combinatorics -> algebra test  ==========" << std::endl;
	std::cout << "(note: expect large values for f and g)" << std::endl;
	Geex::test_algebra(filename_mesh, filename_pts, p);

	return 0;
}