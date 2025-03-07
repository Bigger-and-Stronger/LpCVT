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

// int main(int argc, char** argv) {
// 	std::string filename_mesh = "../data/three_holes.obj";
// 	std::string filename_pts = "../data/three_holes.pts";
// 	std::string output_path;
// 	bool volume = false;
// 	unsigned int p = 2;
// 	bool optimize_mode = false;
//
// 	CLI::App app{"LpCVT"};
// 	argv = app.ensure_utf8(argv);
//
// 	app.add_option("-p,-P", p, "Set L_p, p must be even");
// 	app.add_flag("-v,-V", volume, "Volume or not?");
// 	app.add_flag("-o,-O", optimize_mode, "Optimize?");
// 	app.add_option("meshPath", filename_mesh, "Reference mesh path")->check(CLI::ExistingFile)->required();
// 	app.add_option("ptsPath", filename_pts, "Pts path")->check(CLI::ExistingFile)->required();
// 	app.add_option("outputDir", output_path, "Output dir");
//
// 	CLI11_PARSE(app, argc, argv);
//
// 	/** Process **/
// 	if (p % 2 != 0) {
// 		std::cout << "error: p must be even" << std::endl;
// 		return 1;
// 	}
// 	if (optimize_mode) {
// 		Geex::Mesh M;
// 		M.load(filename_mesh);
//
// 		std::vector<Geex::vec3> pts;
// 		if (!Geex::load_pts(filename_pts, pts))
// 			return 1;
//
// 		Geex::optimize(LBFGSppSolver::LBFGSpp_solver(&M, pts, volume, p));
//
// 		Geex::write_pts(output_path + "res.pts", pts);
// 		filename_pts = output_path + "res.pts";
// 	}
//
// 	std::cout << "============= geometry -> combinatorics test ==========" << std::endl;
// 	Geex::test_combinatorics(filename_mesh, filename_pts, output_path);
// 	std::cout << std::endl;
// 	std::cout << "============= combinatorics -> algebra test  ==========" << std::endl;
// 	std::cout << "(note: expect large values for f and g)" << std::endl;
// 	Geex::test_algebra(filename_mesh, filename_pts, p);
//
// 	return 0;
// }

int main() {
	const std::string input_mesh_path = "../data/three_holes.obj";
	const std::string input_pts_path = "../data/three_holes_surf_rand.pts";
	const std::string output_pts_path = "../data/res_sample_pts.obj";

	// generate random sample pts
	if constexpr (false) {
		Geex::Mesh mesh;
		mesh.load(input_mesh_path);//input_mesh_path);
		std::vector<Geex::vec3> pts;
		Geex::sample_mesh_pts(mesh, pts, 50000);
		Geex::write_pts("../data/test_smp.obj", pts);

		return 1;
	}
	else if constexpr (true) {
		/** load mesh **/
		Geex::Mesh mesh;
		mesh.load(input_mesh_path);

		/** load sample points **/
		std::vector<Geex::vec3> pts;
		load_pts(input_pts_path, pts);

		/** setting params **/
		constexpr bool volume = false;
		const unsigned int p = 8;

		/** optimize **/
		Geex::optimize(LBFGSppSolver::LBFGSpp_solver(&mesh, pts, volume, p));

		/** output **/
		Geex::write_pts(output_pts_path, pts);
	}

	std::cout << std::endl;
	VERBOSE("============= geometry -> combinatorics test ==========");
	Geex::test_combinatorics(input_mesh_path, output_pts_path, "../data/");
	/*
	std::cout << std::endl;
	VERBOSE("============= combinatorics -> algebra test  ==========");
	VERBOSE("(note: expect large values for f and g)");
	Geex::test_algebra(input_mesh_path, output_pts_path);*/
	return 1;
}

/*
int main(int argc, char** argv) {
	if (argc != 3) {
		std::cerr << "usage  : " << argv[0] << " mesh_filename pts_filename" << std::endl;
		std::cerr << "example: " << argv[0] << " data/three_holes.obj data/three_holes.pts" << std::endl;
		return -1;
	}
	std::cerr << "============= geometry->combinatorics test ==========" << std::endl;
	Geex::test_combinatorics(argv[1], argv[2]);
	std::cerr << "============= combinatorics->algebra test  ==========" << std::endl;
	std::cerr << "(note: expect large values for f and g)" << std::endl;
	Geex::test_algebra(argv[1], argv[2]);
	return 0;
}*/
