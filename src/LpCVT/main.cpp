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

#include <LpCVT/combinatorics/delaunay.h>

int main() {
	const std::string input_mesh_path = "../data/three_holes.obj";
	const std::string input_pts_path = "../data/test_sample_pts.obj";
	const std::string output_pts_path = "../data/res_sample_pts.obj";

	// generate random sample pts
	if constexpr (false) {
		Geex::Mesh mesh;
		mesh.load(input_mesh_path);
		std::vector<Geex::vec3> pts;
		Geex::sample_mesh_pts(mesh, pts, 5000);
		Geex::write_pts("../data/test_sample_pts.obj", pts);

		return 1;
	}
	else if constexpr (true) {
		Geex::Optimizer OP;
		OP.read_mesh(input_mesh_path);
		if (!OP.read_pts(input_pts_path))
			return 0;
		OP.optimize();
		OP.write_pts(output_pts_path);
	}

	std::cout << std::endl;
	VERBOSE("============= geometry -> combinatorics test ==========");
	Geex::test_combinatorics(input_mesh_path, output_pts_path, "../data/");
	std::cout << std::endl;
	VERBOSE("============= combinatorics -> algebra test  ==========");
	VERBOSE("(note: expect large values for f and g)");
	Geex::test_algebra(input_mesh_path, output_pts_path);
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
