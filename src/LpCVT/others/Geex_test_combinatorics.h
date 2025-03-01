//
// This file is derived from LpCVT/main.cpp in the original code from https://xueyuhanlang.github.io
// Created by Canjia Huang <canjia7@gmail.com> on 25-3-1.
//

//===================================================================================================
//
// Geometry and combinatorics tests
//
//  Computes and output to a file:
//     Restricted Voronoi Diagram (intersection between a mesh and a 3D Voronoi diagram)
//     Restricted Delaunay Triangulation (dual of a Restricted Voronoi Diagram)
//     Clipped Voronoi Diagram (intersection between the interior of a mesh and a 3D Voronoi diagram)
//
//===================================================================================================

#ifndef GEEX_TEST_COMBINATORICS_H
#define GEEX_TEST_COMBINATORICS_H

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/clipped_VD.h>
#include <LpCVT/algebra/F_Lp.h>
#include <fstream>

namespace Geex {
    /** Loads points from file.
     * file format: "v" << x << y << z **/
    bool load_pts(const std::string& filename, std::vector<vec3>& pts);

    /** Used by save_RDT(). **/
    class SavePrimalTriangle {
    public:
        explicit SavePrimalTriangle(std::ofstream& out): out_(&out) {}

        void operator()(unsigned int i, unsigned int j, unsigned int k) const {
            (*out_) << "f " << i + 1 << " " << j + 1 << " " << k + 1 << std::endl;
        }
    private:
        std::ofstream* out_;
    };

    /** Given a Restricted Voronoi Diagram, saves the Restricted Delaunay
     * Triangulation to a file in alias|wavefront .obj format. **/
    void save_RDT(RestrictedVoronoiDiagram& RVD, const std::string& filename);

    /** used by save_RVD() **/
    class SaveRVDFacets {
    public:
        explicit SaveRVDFacets(std::ostream& out) : out_(out), cur_v_(1), cur_f_(1) {
            out << "# attribute chart facet integer" << std::endl;
        }
        void operator()(unsigned int iv, Mesh* M) const;
    private:
        std::ostream& out_;
        mutable unsigned int cur_v_;
        mutable unsigned int cur_f_;
    };

    /** Saves a Restricted Voronoi Diagram to a file in alias|wavefront .obj format
     * (with Graphite extensions: facet attributes, rename as .eobj to display). **/
    void save_RVD(RestrictedVoronoiDiagram& RVD, const std::string& filename);

    /** Used by save_clippedVD() **/
    class SaveClippedVDFacets {
    public:
        SaveClippedVDFacets(Delaunay* delaunay, std::ostream& out, double shrink) : delaunay_(delaunay), out_(out), shrink_(shrink), cur_(1) {
            out << "# attribute chart facet integer" << std::endl;
        }
        void operator()(unsigned int i, int j, const vec3& p1, const vec3& p2, const vec3& p3) const;
    private:
        Delaunay* delaunay_;
        std::ostream& out_;
        double shrink_;
        mutable unsigned int cur_;
    };

    /** Saves a Clipped Voronoi Diagram to a file in alias|wavefront .obj format
     * (with Graphite extensions: facet attributes, rename as .eobj to display).
     * The 'shrink' parameter is to generate the clipped Voronoi cells as in
     * Figures 4 and 5 in the paper. */
    void save_clippedVD(ClippedVoronoiDiagram& CVD, const std::string& filename, double shrink);

    /** Just used for executable file **/
    void test_combinatorics(const std::string& mesh_filename, const std::string& pts_filename,
        const std::string& output_path);
}

#endif //GEEX_TEST_COMBINATORICS_H
