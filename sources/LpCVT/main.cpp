/*
 * Information removed for anonymous submission
 */ 

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/clipped_VD.h>
#include <LpCVT/algebra/F_Lp.h>
#include <LpCVT/common/line_stream.h>
#include <fstream>

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

namespace Geex {

    //==================================================================================

    /**
     * Loads points from file.
     */

    void load_pts(const std::string& filename, std::vector<vec3>& pts) {
        pts.clear() ;
        std::ifstream in_stream(filename.c_str()) ;
        if(!in_stream) { 
            std::cerr << "Could not open " << filename << std::endl ;
            return ; 
        }
        LineInputStream in(in_stream) ;
        while(!in.eof()) {
            in.get_line() ;
            std::string kw ;
            in >> kw ;
            if(kw == "v") {
                vec3 v ;
                in >> v ;
                pts.push_back(v) ;
            }
        }
    }


    //==================================================================================

    /**
     * Used by save_RDT().
     */
    class SavePrimalTriangle {
    public:
        SavePrimalTriangle(
            std::ofstream& out
        ) : out_(&out) { 
        }
        void operator()(unsigned int i, unsigned int j, unsigned int k) const {
            (*out_) << "f " << i+1 << " " << j+1 << " " << k+1 << std::endl ;
        }

    private:
        std::ofstream* out_ ;
    } ;

    /**
     * Given a Restricted Voronoi Diagram, saves the Restricted Delaunay
     * Triangulation to a file in alias|wavefront .obj format.
     */
    void save_RDT(
        RestrictedVoronoiDiagram& RVD, const std::string& filename
    ) {
        std::cerr << "Computing and saving RDT to " << filename << std::endl ;
        std::ofstream out(filename.c_str()) ;
        for(unsigned int i=0; i<RVD.delaunay()->nb_vertices(); i++) {
            out << "v " << RVD.delaunay()->vertex(i) << std::endl ;
        }
        RVD.for_each_primal_triangle(SavePrimalTriangle(out)) ;
        out.close();
        std::cerr << "Done." << std::endl ;
    }

    //==================================================================================

    /**
     * used by save_RVD()
     */
    class SaveRVDFacets {
    public:
        SaveRVDFacets(
            std::ostream& out
        ) : out_(out), cur_v_(1), cur_f_(1) {
            out << "# attribute chart facet integer" << std::endl ;
        }
        void operator()(unsigned int iv, Mesh* M) const {
            for(unsigned int f=0; f<M->nb_facets(); f++) {
                for(unsigned int i=M->facet_begin(f); i<M->facet_end(f); i++) {
                    const vec3& v = M->vertex(i) ;
                    out_ << "v " << v << std::endl ;
                }
                out_ << "f " ;
                for(unsigned int i=M->facet_begin(f); i<M->facet_end(f); i++) {
                    out_ << cur_v_ << " ";
                    cur_v_++ ;
                }
                out_ << std::endl ;
                out_ << "# attrs f " << cur_f_ << " " << iv << std::endl ;
                cur_f_++ ;
            }
        }        
    private:
        mutable std::ostream& out_ ;
        mutable unsigned int cur_v_ ;
        mutable unsigned int cur_f_ ;
    } ;

    /**
     * Saves a Restricted Voronoi Diagram to a file in alias|wavefront .obj format
     * (with Graphite extensions: facet attributes, rename as .eobj to display).
     */
    void save_RVD(RestrictedVoronoiDiagram& RVD, const std::string& filename) {
        std::ofstream out(filename.c_str()) ;
        if(!out) {
            std::cerr << "could not open file." << std::endl ;
            return ;
        }
        std::cerr << "Computing and saving RVD" << std::endl ;
        bool sym_backup = RVD.symbolic() ;
        RVD.set_symbolic(true) ;
        RVD.for_each_facet(SaveRVDFacets(out)) ;
        RVD.set_symbolic(sym_backup) ;
        std::cerr << "Saved RVD in " << filename << std::endl ;
    }

    //==================================================================================

    /**
     * Used by save_clippedVD()
     */
    class SaveClippedVDFacets {
    public:
        SaveClippedVDFacets(
            Delaunay* delaunay, std::ostream& out, double shrink
        ) : delaunay_(delaunay), out_(out), shrink_(shrink), cur_(1) {
            out << "# attribute chart facet integer" << std::endl ;
        }
        void operator()(unsigned int i, int j, const vec3& p1, const vec3& p2, const vec3& p3) const {
            vec3 x0 = delaunay_->vertex(i) ;
            out_ << "v " << x0 + shrink_ * (p1 - x0) << std::endl ;
            out_ << "v " << x0 + shrink_ * (p2 - x0) << std::endl ;
            out_ << "v " << x0 + shrink_ * (p3 - x0) << std::endl ;
            out_ << "f " << cur_ << " " << cur_ + 1 << " " << cur_ + 2 << std::endl ;
            cur_ += 3 ;
            out_ << "# attrs f " << ((cur_ - 1) / 3) << " " << i << std::endl ;
        }
    private:
        Delaunay* delaunay_ ;
        mutable std::ostream& out_ ;
        double shrink_ ;
        mutable unsigned int cur_ ;
    } ;

    /**
     * Saves a Clipped Voronoi Diagram to a file in alias|wavefront .obj format
     * (with Graphite extensions: facet attributes, rename as .eobj to display).
     * The 'shrink' parameter is to generate the clipped Voronoi cells as in
     * Figures 4 and 5 in the paper.
     */
    void save_clippedVD(ClippedVoronoiDiagram& CVD, const std::string& filename, double shrink) {
        std::ofstream out(filename.c_str()) ;
        if(!out) {
            std::cerr << "could not open file." << std::endl ;
            return ;
        }
        std::cerr << "Computing and saving clipped VD" << std::endl ;
        CVD.for_each_triangle(SaveClippedVDFacets(CVD.delaunay(), out,shrink)) ;
        std::cerr << "Saved clipped VD in " << filename << std::endl ;
    }


    //==================================================================================

    void test_combinatorics(const std::string& mesh_filename, const std::string& pts_filename) {
        Mesh M ;
        unsigned int nb_borders = M.load(mesh_filename) ;
        std::vector<vec3> pts ;
        load_pts(pts_filename, pts) ;
        Delaunay* delaunay = Delaunay::create("CGAL") ;
        RestrictedVoronoiDiagram RVD(delaunay, &M) ;
        ClippedVoronoiDiagram CVD(delaunay, &M) ;

        delaunay->set_vertices(pts) ;
        save_RVD(RVD, "rvd.obj") ;
        save_RDT(RVD, "rdt.obj") ;
        if(nb_borders == 0) {
           save_clippedVD(CVD, "cvd.obj", 0.7) ;
        }
        delete delaunay ;
    }

}

//===================================================================================================
//
// Algebra tests
//
//  Computes :
//     F_{L_p} in the surface case
//     F_{L_p} in the volume case
//
//=================================================================================================== 

namespace Geex {

    /**
     * Used by get_combinatorics() in volume mode
     */
    class MemorizeIndices {
    public:
        MemorizeIndices(
            std::vector<int>& I_in,
            std::vector<vec3>& C_in
        ) : I(I_in), C(C_in) {
            I.resize(0) ;
            C.resize(0) ;
        }

        void operator() (
            unsigned int i, 
            int j,
            const VertexEdge& v1, 
            const VertexEdge& v2, 
            const VertexEdge& v3
        ) const {
            I.push_back(i) ;
            I.push_back(v1.sym[2]) ;
            I.push_back(v1.sym[1]) ;
            I.push_back(v1.sym[0]) ;
            I.push_back(v2.sym[2]) ;
            I.push_back(v2.sym[1]) ;
            I.push_back(v2.sym[0]) ;
            I.push_back(v3.sym[2]) ;
            I.push_back(v3.sym[1]) ;
            I.push_back(v3.sym[0]) ;
            C.push_back(v1) ;
            C.push_back(v2) ;
            C.push_back(v3) ;
        }
    private:
        mutable std::vector<int>& I ;
        mutable std::vector<vec3>& C ;
    } ;

    /**
     * Used by get_combinatorics() in surface mode
     */
    class MemorizeIndicesAndFacets{
    public:
        MemorizeIndicesAndFacets(
            const RestrictedVoronoiDiagram& RVD_in,
            std::vector<int>& I_in,
            std::vector<vec3>& C_in,
            std::vector<int>& F_in
        ) : RVD(RVD_in), I(I_in), C(C_in), F(F_in) {
            I.resize(0) ;
            C.resize(0) ;
            F.resize(0) ;
        }

        void operator() (
            unsigned int i, 
            const VertexEdge& v1, 
            const VertexEdge& v2, 
            const VertexEdge& v3
        ) const {
            I.push_back(i) ;
            I.push_back(v1.sym[2]) ;
            I.push_back(v1.sym[1]) ;
            I.push_back(v1.sym[0]) ;
            I.push_back(v2.sym[2]) ;
            I.push_back(v2.sym[1]) ;
            I.push_back(v2.sym[0]) ;
            I.push_back(v3.sym[2]) ;
            I.push_back(v3.sym[1]) ;
            I.push_back(v3.sym[0]) ;
            F.push_back(RVD.current_facet()) ;
            C.push_back(v1) ;
            C.push_back(v2) ;
            C.push_back(v3) ;
        }
    private:
        const RestrictedVoronoiDiagram& RVD ;
        mutable std::vector<int>& I ;
        mutable std::vector<vec3>& C ;
        mutable std::vector<int>& F ;
    } ;


    /**
     * Gets the combinatorics of the integration simplices,
     * i.e. 10 integers per integration simplex.
     * (see Section 3.1 in the paper)
     * Returns also the array of C vertices (three per integration simplex).
     * Since they are easy to get during the combinatorial phase, they are
     * computed here and kept for the algebraic phase.
     *
     * In 2D mode (volume = false), returns also the array F.
     *   F[i] indicates the facet that contains the i-th integration simplex.
     *
     */
    void get_combinatorics(
        Mesh* M, const std::vector<vec3>& pts, 
        std::vector<int>& I, std::vector<vec3>& C, std::vector<int>& F, bool volume
    ) {
        Delaunay* delaunay = Delaunay::create("CGAL") ;
        delaunay->set_vertices(pts) ;
        if(volume) {
            ClippedVoronoiDiagram CVD(delaunay,M) ;
            CVD.for_each_triangle(MemorizeIndices(I,C)) ;
        } else {
            RestrictedVoronoiDiagram RVD(delaunay,M) ;
            RVD.set_symbolic(true) ;
            RVD.for_each_triangle(MemorizeIndicesAndFacets(RVD,I,C,F)) ;
        }
        delete delaunay ;
    }

    /**
     * Computes F_{L_p} and its gradient.
     */
    void compute_F_g(Mesh* m, const std::vector<vec3>& pts, unsigned int p, bool volume) {
        std::cerr << "nb pts = " << pts.size() << "   nb facets = " << m->nb_facets() << std::endl ;
        std::vector<int> I ;
        std::vector<vec3> C ;
        std::vector<int> F ;
        get_combinatorics(m, pts, I, C, F, volume) ;
        unsigned int nb_integration_simplices = (unsigned int)I.size() / 10 ;
        std::vector<mat3> M(nb_integration_simplices) ;
        for(unsigned int i=0; i<M.size(); i++) {
            M[i].load_identity() ; 
                // or replace with anisotropy field
                //   In 2D: use F[i] to retreive the index of the facet that contains
                //      the current integration simplex (and access an array of per-facet anisotropy).
                //   In 3D: use geometric search from the centroid of the current
                //      integration simplex.
        }
        std::vector<plane3> Q(m->nb_facets()) ;
        for(unsigned int i=0; i<m->nb_facets(); i++) {
            Q[i] = m->facet_plane(i) ;
        }
        std::vector<double> g(pts.size() * 3) ;
        double f = compute_F_Lp(volume, p, m, I, C, pts, Q, M, g) ;
        double gnorm = 0.0 ;
        for(unsigned int i=0; i<g.size(); i++) {
            gnorm += g[i]*g[i] ;
        }
        gnorm = ::sqrt(gnorm) ;
        std::cerr.precision(16);
        std::cerr << (volume ? "volume " : "surface ") 
                  << "F_L" << p << ":" 
                  << "f=" << std::scientific << f << "  g=" << gnorm << std::endl ;
    }

    void test_algebra(const std::string& mesh_filename, const std::string& pts_filename) {
        Mesh M ;
        unsigned int nb_borders = M.load(mesh_filename) ;
        std::vector<vec3> pts ;
        load_pts(pts_filename, pts) ;
        if(nb_borders == 0) {
	   std::cerr << "          ========== volume LpCVT test ======" << std::endl ;
	   compute_F_g(&M, pts, 4, true) ;
	}
        std::cerr << "          ========== surface LpCVT test ======" << std::endl ;
        compute_F_g(&M, pts, 4, false) ;
    }

}

int main(int argc, char** argv) {
    if(argc != 3) {
        std::cerr << "usage  : " << argv[0] << " mesh_filename pts_filename" << std::endl ;
        std::cerr << "example: " << argv[0] << " data/three_holes.obj data/three_holes.pts" << std::endl ;
        return -1 ;
    }
    std::cerr << "============= geometry->combinatorics test ==========" << std::endl ;
    Geex::test_combinatorics(argv[1], argv[2]) ;
    std::cerr << "============= combinatorics->algebra test  ==========" << std::endl ;
    std::cerr << "(note: expect large values for f and g)" << std::endl ;
    Geex::test_algebra(argv[1], argv[2]) ;
    return 0 ;
}
