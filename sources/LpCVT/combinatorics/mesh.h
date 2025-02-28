/*
 *  Copyright (c) 2010, INRIA, Project ALICE
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://alice.loria.fr
 *
 *     ALICE Project
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */


#ifndef __GEEX_CVT_TOPO_POLY_MESH__
#define __GEEX_CVT_TOPO_POLY_MESH__

#include <LpCVT/common/types.h>
#include <LpCVT/combinatorics/symbolic_vertex.h>
#include <vector>
#include <iostream>

namespace Geex {

    typedef unsigned int TopoFlags ;
    typedef unsigned int TopoFlag ;

    /**
     * Vertex of a polygon + information
     * attached to the edge starting from
     * this vertex.
     */
    class VertexEdge : public vec3Sym {
    public:
        enum {ORIGINAL = 1, INTERSECT = 2} ;

        VertexEdge(
            const VertexEdge& rhs
        ) : vec3Sym(rhs, rhs.sym), flags(rhs.flags), f(rhs.f), w(rhs.w) {
        }

        VertexEdge(
            const vec3Sym& rhs
        ) : vec3Sym(rhs), flags(0), f(-1), w(1.0) {
        }

        VertexEdge(
            const vec3& p, double w_in = 1.0
        ) : vec3Sym(p), flags(0), f(-1), w(w_in) {
        }

        VertexEdge() { }

        void set_point(const vec3& p) {
            vec3::operator=(p) ;
        }

        void set_edge(const VertexEdge& rhs) {
            flags = rhs.flags ;
            f = rhs.f ;
        }

        VertexEdge& operator=(const VertexEdge& rhs) {
            vec3::operator=(rhs) ;
            set_edge(rhs) ;
            w = rhs.w ;
            sym = rhs.sym ;
            return *this ;
        }

        void clear() { flags = 0 ; f = -1 ; w = 1.0 ; }
        void set_flag(TopoFlag f) { flags |= f ; }
        void unset_flag(TopoFlag f) { flags &= ~f ; }
        bool check_flag(TopoFlag f) const { return ((flags & f) != 0) ; }

        TopoFlags flags ;    // Indicates the type of edge 
                             // (virtual, original or intersection).
        int f ;              // The index of the other facet incident 
                             //  to this edge (or -1 on boundary).
        double w ;           // Vertex weight (e.g., curvature). 
    } ;

    /**
     * Information attached to a facet. 
     * flags are used to mark/unmark facets during traversal.
     */
    class FacetInfo {
    public:
        FacetInfo() : flags(0), id(~0) { }
        TopoFlags flags ;
        unsigned int id ;
    } ;

    /**
     * Mesh stores facet connectivity of a polygon mesh. 
     * Vertices are duplicated. Facets store pointers to adjacent
     * facets.
     * Facets are stored using the CRS (compressed row storage) 
     * format, i.e. the vertices of the f-th facet are traversed
     * using:
     * for(unsigned int i=M.facet_begin(f); i<M.facet_end(f); i++) {
     *    do something with M.vertex(i)
     * }
     *
     * In addition, Mesh can be used as a polygon stack
     * (this is used by the RVD computation algorithm).
     */
    class Mesh {
    public:
        Mesh(
        ) : in_facet_(false){
            facet_ptr_.push_back(0) ; 
        }

        Mesh(
            const Mesh& rhs
        ) : in_facet_(false){
            facet_ptr_.push_back(0) ; 
        }

        void clear() {
            // Note: we use resize(0) instead of clear()
            // since this is guaranteed to keep the reserved
            // memory (this is what we want for our "ping-pong"
            // buffers).
            vertex_.resize(0) ; 
            facet_ptr_.resize(0) ; 
            facet_info_.resize(0) ; 
            facet_ptr_.push_back(0) ; 
            original_vertices_.resize(0) ;
            vertex_index_.resize(0) ;
        }

        //==----------------------- Construction

        void begin_facet() { 
            in_facet_ = true ; 
        }

        void add_vertex(const VertexEdge& ve) {
            vertex_.push_back(ve) ;
        }

        void end_facet() {
            facet_ptr_.push_back((unsigned int)vertex_.size()) ;
            facet_info_.push_back(FacetInfo()) ;
            in_facet_ = false ;
        }

        bool in_facet() const { return in_facet_ ; }

        //==----------------------- Access

        unsigned int nb_facets() const { return (unsigned int)facet_ptr_.size() - 1 ; }
        unsigned int nb_vertices() const { return (unsigned int)vertex_.size() ; }
        unsigned int facet_begin(unsigned int f) const {
            return facet_ptr_[f] ;
        }

        unsigned int facet_end(unsigned int f) const {
            return facet_ptr_[f+1] ;
        }

        unsigned int facet_size(unsigned int f) const {
            return facet_end(f) - facet_begin(f) ;
        }

        VertexEdge& facet_vertex(unsigned int f, unsigned int i) {
            return vertex(facet_begin(f)+i) ;
        }

        const VertexEdge& facet_vertex(unsigned int f, unsigned int i) const {
            return vertex(facet_begin(f)+i) ;
        }

        unsigned int next_around_facet(unsigned int f, unsigned int i) const {
            return (i+1 == facet_end(f) ? facet_begin(f) : i+1) ;
        }

        unsigned int prev_around_facet(unsigned int f, unsigned int i) const {
            return (i == facet_begin(f) ? facet_end(f)-1 : i-1) ;
        }

        /**
         * do_it is supposed to overload operator()(
         *   const VertexEdge& v1, 
         *   const VertexEdge& v2, 
         *   const VertexEdge& v3 )
         * Can be for instance an inline function.
         */
        template<class T> void for_each_triangle(unsigned int f, T& do_it) const {
            unsigned int i0 = facet_begin(f) ;
            for(unsigned int i = facet_begin(f)+1; i+1<facet_end(f); i++) {
                do_it(vertex(i0), vertex(i), vertex(i+1)) ;
            }
        }

        /**
         * do_it is supposed to overload operator()(
         *   const VertexEdge& v1, 
         *   const VertexEdge& v2, 
         *   const VertexEdge& v3 )
         * Can be for instance an inline function.
         */
        template<class T> void for_each_triangle(
            unsigned int f, const T& do_it
        ) const {
            unsigned int i0 = facet_begin(f) ;
            for(unsigned int i = facet_begin(f)+1; i+1<facet_end(f); i++) {
                do_it(vertex(i0), vertex(i), vertex(i+1)) ;
            }
        }

        /**
         * do_it is supposed to overload operator()(
         *   const VertexEdge& v1, 
         *   const VertexEdge& v2, 
         *   const VertexEdge& v3 )
         * Can be for instance an inline function.
         */
        template<class T> void for_each_triangle(T& do_it) const {
            for(unsigned int f=0; f<nb_facets(); f++) {
                unsigned int i0 = facet_begin(f) ;
                for(unsigned int i = facet_begin(f)+1; i+1<facet_end(f); i++) {
                    do_it(vertex(i0), vertex(i), vertex(i+1)) ;
                }
            }
        }

        /**
         * do_it is supposed to overload operator()(
         *   const VertexEdge& v1, 
         *   const VertexEdge& v2, 
         *   const VertexEdge& v3 )
         * Can be for instance an inline function.
         */
        template<class T> void for_each_triangle(const T& do_it) const {
            for(unsigned int f=0; f<nb_facets(); f++) {
                unsigned int i0 = facet_begin(f) ;
                for(unsigned int i = facet_begin(f)+1; i+1<facet_end(f); i++) {
                    do_it(vertex(i0), vertex(i), vertex(i+1)) ;
                }
            }
        }

        // Note: this is not the centroid, but this is 
        // OK for what we do with it.
        vec3 facet_center(unsigned int f) const {
            vec3 result(0.0, 0.0, 0.0) ;
            double n(0.0) ;
            for(unsigned int i=facet_begin(f); i<facet_end(f); i++) {
                result += vertex(i) ;
                n++ ;
            }
            return (1.0/n)*result ;
        }

        double facet_area(unsigned int f) const {
            double result = 0.0 ;
            unsigned int i0 = facet_begin(f) ;
            for(unsigned int i = facet_begin(f)+1; i+1<facet_end(f); i++) {
                vec3 v1 = vertex(i) - vertex(i0) ;
                vec3 v2 = vertex(i+1) - vertex(i0) ;
                result += cross(v1,v2).length() ;
            }
            result /= 2.0 ;
            return result ;
        }
        
        vec3 facet_normal(unsigned int f) const {
            vec3 result = cross(
                facet_vertex(f,1) - facet_vertex(f,0), facet_vertex(f,2) - facet_vertex(f,0)
            ) ; 
            int d = facet_size(f) ;
            // Check for degeneracies
            if(d > 3) {
                for(unsigned int i=1; i<facet_size(f); i++) {
                    int i1 = i ;
                    int i2 = (i1+1)%d ;
                    int i3 = (i2+1)%d ;
                    result += cross(
                        facet_vertex(f,i2) - facet_vertex(f,i1), facet_vertex(f,i3) - facet_vertex(f,i1)
                    ) ;                    
                }
            }
            return result ;
        }

        plane3 facet_plane(unsigned int f) const { 
            return plane3(facet_vertex(f,0), facet_normal(f)) ; 
        }

        VertexEdge& vertex(unsigned int v) {
            return vertex_[v] ;
        }

        const VertexEdge& vertex(unsigned int v) const {
            return vertex_[v] ;
        }

        //==----------------------- Marking

        void mark_facet(unsigned int f) { 
            facet_info_[f].flags = 1 ;
        }
        void unmark_facet(unsigned int f) { 
            facet_info_[f].flags = 0 ;
        }
        bool facet_is_marked(unsigned int f) const { 
            return (facet_info_[f].flags != 0) ;
        }
        const FacetInfo& facet_info(unsigned int f) const {
            return facet_info_[f] ;
        }
        FacetInfo& facet_info(unsigned int f) {
            return facet_info_[f] ;
        }

        //==----------------------- Stack behavior

        unsigned int top() const { 
            return nb_facets() - 1 ; 
        }
        
        VertexEdge& top_vertex() {
            return *(vertex_.rbegin()) ;
        }

        const VertexEdge& top_vertex() const {
            return *(vertex_.rbegin()) ;
        }

        bool empty() { return nb_facets() == 0 ; }

        // copies facet to target.
        void copy_facet(unsigned int f, Mesh& target) {
            target.begin_facet() ;
            for(unsigned int v = facet_begin(f) ; v < facet_end(f) ; v++) {
                target.add_vertex(vertex(v)) ;
            }
            target.end_facet() ;
            target.facet_info(target.top()) = facet_info(f) ;
        }

        void copy_top_facet(Mesh& target) {
            copy_facet(top(), target) ;
        }

        // removes last facet
        void pop_facet() {
            facet_ptr_.pop_back() ;
            facet_info_.pop_back() ;
            unsigned int nb_v = *(facet_ptr_.rbegin()) ;
            while(vertex_.size() > nb_v) {
                vertex_.pop_back() ;
            }
        }


        // copies last facet to target and
        // removes last facet.
        void pop_facet(Mesh& target) {
            copy_top_facet(target) ;
            pop_facet() ;
        }

        //==----------------------- IO

        // Returns the number of borders.
        unsigned int load(const std::string& filename) ;

        // Note: does not save adjacencies, just geometry.
        void save(const std::string& filename) ;

        //==----------------------- Clipping

        void intersect(
            VertexEdge& I, 
            const VertexEdge& v1, const VertexEdge& v2, 
            const plane3& P
        ) const {
            double h = fabs( v1.x * P.a + v1.y * P.b + v1.z * P.c + P.d );
            double l = fabs( v2.x * P.a + v2.y * P.b + v2.z * P.c + P.d );
            double hl = h + l;
            //if hl is very small, hl > 0 will be false
            if (hl > 0) {
                double t_h = h / hl;
                double t_l = l / hl;
                I = t_l * v1 + t_h * v2;
                I.w = t_l * v1.w + t_h * v2.w;
            } else {
                I = 0.5 * (v1 + v2);
                I.w = 0.5 * (v1.w + v2.w);
            }
        }

        void intersect(
            VertexEdge& I, 
            const VertexEdge& v1, const VertexEdge& v2, 
            const plane3& P,
            unsigned int E
        ) const {
            double h = fabs( v1.x * P.a + v1.y * P.b + v1.z * P.c + P.d );
            double l = fabs( v2.x * P.a + v2.y * P.b + v2.z * P.c + P.d );
            double hl = h + l;
            //if hl is very small, hl > 0 will be false
            if (hl > 0) {
                double t_h = h / hl;
                double t_l = l / hl;
                I = t_l * v1 + t_h * v2;
                I.w = t_l * v1.w + t_h * v2.w;
            } else {
                I = 0.5 * (v1 + v2);
                I.w = 0.5 * (v1.w + v2.w);
            }
            sets_intersect(v1.sym, v2.sym, I.sym) ;
            I.sym.add_bisector(E) ;
        }

        // ------------------ Queries ---------------------------------

        
        // Returns in w1 and w2 the two extremities
        // of the edge shared by facets q1 and q2.
        inline void find_edge_extremities(
            unsigned int q1, unsigned int q2,
            vec3& w1, vec3& w2
        ) const {
            if(q1 > q2) { 
                unsigned int tmp = q1 ; q1 = q2; q2 = tmp ;
            }
            const vec3* v1 = 0 ; const vec3* v2 = 0 ;
            if(q2 >= nb_facets()) {
                // q2 is a "Virtual facet" on the boundary, indicates the index of
                // the edge in the other facet q1.
                // Note: since w.sym is *ordered*, this will be always q2 that is
                // virtual and q1 that is a facet index.
                assert(q1 < nb_facets()) ;
                unsigned int iv1 = facet_begin(q1) + q2 - nb_facets() ;
                v1 = &vertex(iv1) ;
                v2 = &vertex(next_around_facet(q1, iv1)) ;
            } else {
                // Find the edge (v1, v2) of facet q1 adjacent to facet q2
                for(unsigned int i=facet_begin(q1); i<facet_end(q1); i++) {
                    if(vertex(i).f == q2) {
                        v1 = &vertex(i) ;
                        v2 = (i == facet_end(q1) - 1) ? &(vertex(facet_begin(q1))) : &(vertex(i+1)) ;
                        break ;
                    }
                }
            }
            assert(v1 != 0 && v2 != 0) ;
            w1 = *v1 ; w2 = *v2 ;
        }

        // Returns the vertex shared by facets q1,q2 and q3
        inline const vec3& find_vertex(
            unsigned int q1, unsigned int q2, unsigned int q3
        ) const {

            for(unsigned int i=facet_begin(q1); i<facet_end(q1); i++) {
                unsigned int j =next_around_facet(q1,i) ;
                const VertexEdge& v1 = vertex(i) ;
                const VertexEdge& v2 = vertex(j) ;
                if(
                    v1.f == q2 && v2.f == q3 ||
                    v2.f == q2 && v1.f == q3
                ) {
                    return v2 ;
                }
            }

            for(unsigned int i=facet_begin(q2); i<facet_end(q2); i++) {
                unsigned int j =next_around_facet(q2,i) ;
                const VertexEdge& v1 = vertex(i) ;
                const VertexEdge& v2 = vertex(j) ;
                if(
                    v1.f == q3 && v2.f == q1 ||
                    v2.f == q3 && v1.f == q1
                ) {
                    return v2 ;
                }
            }

            for(unsigned int i=facet_begin(q3); i<facet_end(q3); i++) {
                unsigned int j =next_around_facet(q3,i) ;
                const VertexEdge& v1 = vertex(i) ;
                const VertexEdge& v2 = vertex(j) ;
                if(
                    v1.f == q1 && v2.f == q2 ||
                    v2.f == q1 && v1.f == q2
                ) {
                    return v2 ;
                }
            }
             
            return *(vec3*)0 ;
        }


        // ------------------ METIS interface -------------------------

        // Note: only works if compiled with METIS support (-DWITH_METIS)
        void partition(
            unsigned int nb_parts, std::vector<Mesh>& parts
        ) const ;

        // ------------------ Additional information

        const std::vector<vec3>& original_vertices() const {
            return original_vertices_;
        }

        //   Vertices are duplicated in the facets,
        // but they have a single vertex index (that
        // is valid for meshes loaded from a file).
        unsigned int vertex_index(unsigned int i) const {
            return vertex_index_[i] ;
        }

        // ---------- Signed volume and orientation

        double area() const ; 
        double signed_volume() const ;
        bool orientation() const { return orientation_ ; }

    protected:
        void init_symbolic_vertices() ;

    private:
        bool in_facet_ ;
        std::vector<VertexEdge> vertex_ ;
        std::vector<unsigned int> facet_ptr_ ;
        std::vector<FacetInfo> facet_info_ ;
        std::vector<vec3> original_vertices_ ;
        std::vector<unsigned int> vertex_index_ ;
        bool orientation_ ;

    private:
        Mesh& operator=(const Mesh& rhs) ;
    } ;

}

#endif
