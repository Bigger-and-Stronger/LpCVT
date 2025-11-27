
#ifndef __GEEX_CVT_RVD_PREDICATES__
#define __GEEX_CVT_RVD_PREDICATES__

#include <LpCVT/common/types.h>

namespace Geex {

    class VertexEdge ;
    class Mesh ;
    class RestrictedVoronoiDiagram ;
    class vec3Sym ;

    class RVD_predicates {
    public:
        /*
         * returns the side of vertex v with respect
         * to the bisector of [p1,p2].
         * Positive side (resp. negative) is p1's side (resp. p2)
         * v needs to have symbolic information (RVD needs
         * to be in symbolic or exact mode).
         */
        static Sign side(
            const vec3& p1, const vec3& p2,
            const vec3Sym& v,
            const RestrictedVoronoiDiagram* RVD
        ) ;

        static void set_verbose(bool x) ;
        static void begin_stats() ;
        static void end_stats() ;

        /**
	 * Given the symbolic representation of a vertex, returns
	 * its (inexact) geometry. This function is used to debug 
	 * code that manipulates symbolic information.
	 */ 
        static vec3 sym_to_inexact_geometry(
            unsigned int center_vertex, const vec3Sym& v,
            const RestrictedVoronoiDiagram* RVD
        ) ;

    protected:
        static bool find_edge(
            const Mesh* M, int f1, int f2, vec3& q1, vec3& q2
        ) ;

    } ;

}


#endif

