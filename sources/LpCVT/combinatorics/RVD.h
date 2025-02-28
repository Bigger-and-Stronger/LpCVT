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

#ifndef __GEEX_CVT_RVD__
#define __GEEX_CVT_RVD__

#include <LpCVT/combinatorics/mesh.h>
#include <LpCVT/combinatorics/delaunay_skel.h>
#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/exact/RVD_predicates.h>
#include <stack>
#include <set>

namespace Geex {

    /*
     * Returns a seed near a given facet. Does not
     * need to be exactly the nearest seed.
     */
    inline unsigned int find_seed_near_facet(
        const Delaunay* delaunay, const Mesh* M, unsigned int f
    ) {
        vec3 p = M->facet_center(f) ;
        return delaunay->nearest_vertex_id(p) ;
    }

    //--------------------------------------------------------------------------

    // Used by RVD internally to implement for_each_triangle()
    template <class T> class TriangleAction {
    public:
        TriangleAction(const T& do_it_in) : do_it(do_it_in) { }
        void operator()(unsigned int v, Mesh* M) const {
            for(unsigned int f=0; f<M->nb_facets(); f++) {
                unsigned int i0 = M->facet_begin(f) ;
                for(unsigned int i = M->facet_begin(f)+1; i+1<M->facet_end(f); i++) {
                    do_it(v, M->vertex(i0), M->vertex(i), M->vertex(i+1)) ;
                }
            }
        }
    private:
        const T& do_it ;
    } ;

    // Used by RVD internally to implement for_each_halfedge()
    template <class T> class HalfedgeAction {
    public:
        HalfedgeAction(const T& do_it_in) : do_it(do_it_in) { }
        void operator()(unsigned int v, Mesh* M) const {
            for(unsigned int f=0; f<M->nb_facets(); f++) {
                for(unsigned int i = M->facet_begin(f); i<M->facet_end(f); i++) {
                    if(M->vertex(i).check_flag(VertexEdge::INTERSECT)) {
                        do_it(v, M->vertex(i), M->vertex(M->next_around_facet(f,i))) ;
                    }
                }
            }
        }
    private:
        const T& do_it ;
    } ;

    // Used by RVD internally to implement for_each_primal_triangle()
    template <class T> class PrimalTriangleAction {
    public:
        PrimalTriangleAction(const T& do_it_in) : do_it(do_it_in) { }
        void operator()(unsigned int iv1, Mesh* M) const {
            for(unsigned int f=0; f<M->nb_facets(); f++) {
                for(unsigned int i = M->facet_begin(f); i<M->facet_end(f); i++) {
                    const VertexEdge& ve = M->vertex(i) ;
                    // Primal triangles correspond to vertices of
                    // the RVD that are on two bisectors.  
                    unsigned int nb_bisect = ve.sym.nb_bisectors() ;
                    if(nb_bisect >= 2) {
                        bool found = false ;
                        for(unsigned int k1=0; k1<nb_bisect-1 && !found; k1++) {
                            for(unsigned int k2=k1+1; k2<nb_bisect && !found; k2++) {
                                unsigned int iv2 = ve.sym.bisector(k1) ; 
                                unsigned int iv3 = ve.sym.bisector(k2) ;
                                // Make sure each triangle is generated once only
                                // (skip the triangle if iv1 is not the smallest index)
                                if(iv1 < iv2 && iv1 < iv3) {
                                    // Check whether the orientation of the triangle is
                                    // consistent with the orientation of the facet. If
                                    // its not the case, swap iv2 and iv3
                                    unsigned int j1 = M->prev_around_facet(f,i) ;
                                    unsigned int j2 = M->next_around_facet(f,i) ;

                                    const VertexEdge& ve2 = M->vertex(j1) ; 
                                    const VertexEdge& ve3 = M->vertex(j2) ;

                                    if(ve2.sym.has_bisector(iv2) && ve3.sym.has_bisector(iv3)) {
                                        do_it(iv1, iv2, iv3) ;
                                        found = true ;
                                    } else if(ve2.sym.has_bisector(iv3) && ve3.sym.has_bisector(iv2)) {
                                        do_it(iv1, iv3, iv2) ;
                                        found = true ;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    private:
        const T& do_it ;
    } ;

    //--------------------------------------------------------------------------


    struct RVDStackItem {
        RVDStackItem(
            unsigned int f_in, unsigned int cell_in
        ) : f(f_in), cell(cell_in) { }
        unsigned int f ;
        unsigned int cell ;
    } ;

    typedef std::stack< RVDStackItem > RVDStack ;
    typedef std::vector<unsigned int> TopoCellSet ;


    /*
     * Computes the restricted voronoi diagram using propagation.
     * This version uses simple cell clipping.
     */
    class RestrictedVoronoiDiagram {
        typedef RestrictedVoronoiDiagram thisclass ;

    protected:
        void init_visited(Mesh* M) {
            m = M ;
            if(M->nb_facets() != visited_size_) {
                delete[] visited_ ; 
                visited_size_ = M->nb_facets() ;
                visited_ = new TopoCellSet[visited_size_] ;
            } else {
                for(unsigned int i=0; i<visited_size_; i++) {
                    // We use resize(0) instead of clear() since
                    // this preserves the allocated memory.
                    visited_[i].resize(0) ;
                }
            }
        }


    public:
        RestrictedVoronoiDiagram() : 
            visited_(0), visited_size_(0), m(0), delaunay_(0), symbolic_(false), exact_(false){ 
        }
        
        RestrictedVoronoiDiagram(Delaunay* delaunay_in, Mesh* m_in) : 
            visited_(0), visited_size_(0), m(m_in), delaunay_(delaunay_in), 
            symbolic_(false), exact_(false) { 
        }
        
        ~RestrictedVoronoiDiagram() { delete[] visited_ ; visited_ = 0 ; }

        // In symbolic mode, each generated vertex knows the ID's of the planes
        // it belongs to.
        void set_symbolic(bool x) { symbolic_ = x ; }

        bool symbolic() const { return symbolic_ ; }

        // In exact mode, all combinatorial decisions (predicates) use exact
        // arithmetics.
        void set_exact(bool x) { exact_ = x ; }

        bool exact() const { return exact_ ; }

        const Mesh* mesh() const { return m ; }
        Mesh* mesh() { return m ; }

        const Delaunay* delaunay() const { return delaunay_ ; }
        Delaunay* delaunay() { return delaunay_ ; }

        // Low-level API. Client code may use for_each_facet(), for_each_triangle() or
        // for_each_primal_triangle() instead.
        template <class DEL, class ACTION> inline void 
        compute(
            Mesh* M, DEL* delaunay, const DelaunaySkeleton* skel, const ACTION& action
        ) {
            current_mesh_ = 0 ;
            init_visited(M) ;
            for(unsigned int i=0; i<M->nb_facets(); i++) {
                if(!M->facet_is_marked(i)) {
                    unsigned int cell = find_seed_near_facet(delaunay, M, i) ;
                    push(i,cell,false) ;
                    while(!S.empty()) {
                        Mesh* ping = &M1 ;
                        Mesh* pong = &M2 ;
                        unsigned int f = S.top().f ;
                        unsigned int cell = S.top().cell ;

                        current_facet_ = f ;

                        S.pop() ;
                        ping->clear() ;
                        M->copy_facet(f, *ping) ;
                        
                        // Compute triangle clipped by cell and
                        // send neighbors to stack
                        clip_by_cell(skel, cell, ping, pong) ;
                        current_mesh_ = ping ;
                        // Apply action to triangles intersected by cell.
                        action(cell, ping) ;
                        
                        for(unsigned int i=0; i<ping->nb_vertices(); i++) {
                            int neigh_f = ping->vertex(i).f ;
                            if(neigh_f >= 0 && !M->facet_is_marked(neigh_f)) {
                                push(neigh_f, cell) ;
                            }
                        }
                    }
                }
            }

            // Final cleanup: unmark all triangles
            for(unsigned int f=0; f<M->nb_facets(); f++) {
                M->unmark_facet(f) ;
            }
            current_mesh_ = 0 ;
        }

        // ACTION needs to implement 
        //      operator()(unsigned int c, Mesh* M) const
        //   where c denotes the index of the current Voronoi cell (or Delaunay vertex).
        template <class ACTION> inline void for_each_facet(const ACTION& action) {
            if(exact_) {
                RVD_predicates::begin_stats() ;
            }
            this->template compute<Delaunay, ACTION>(m, delaunay_, delaunay_->skeleton(), action) ;
            if(exact_) {
                RVD_predicates::end_stats() ;
            }
        }

        // ACTION needs to implement:
        //      operator()(unsigned int c, const VertexEdge& v1, v2, v3) const
        //   where c denotes the index of the current Voronoi cell (or Delaunay vertex).
        template <class TRIACTION> inline void for_each_triangle(const TRIACTION& action) {
            if(exact_) {
                RVD_predicates::begin_stats() ;
            }
            this->template compute<Delaunay, TriangleAction<TRIACTION> >(
                m, delaunay_, delaunay_->skeleton(), TriangleAction<TRIACTION>(action)
            ) ;
            if(exact_) {
                RVD_predicates::end_stats() ;
            }
        }

        // ACTION needs to implement:
        //      operator()(unsigned int c, const VertexEdge& v1, v2) const
        //   where c denotes the index of the current Voronoi cell (or Delaunay vertex).
        template <class HEACTION> inline void for_each_halfedge(const HEACTION& action) {
            if(exact_) {
                RVD_predicates::begin_stats() ;
            }
            this->template compute<Delaunay, HalfedgeAction<HEACTION> >(
                m, delaunay_, delaunay_->skeleton(), HalfedgeAction<HEACTION>(action)
            ) ;
            if(exact_) {
                RVD_predicates::end_stats() ;
            }
        }


        // ACTION needs to implement:
        //      operator()(unsigned int i, unsigned j, unsigned int k) const
        //   where i,j,k denote the three indices of the Delaunay vertices 
        //   that define the primal triangle.
        template <class PRIMTRIACTION> inline void for_each_primal_triangle(
            const PRIMTRIACTION& action
        ) {
            if(exact_) {
                RVD_predicates::begin_stats() ;
            }
            // PrimalTriangleAction needs symbolic mode
            bool sym_backup = symbolic_ ;
            symbolic_ = true ;
            this->template compute<Delaunay, PrimalTriangleAction<PRIMTRIACTION> >(
                m, delaunay_, delaunay_->skeleton(), PrimalTriangleAction<PRIMTRIACTION>(action)
            ) ;
            symbolic_ = sym_backup ;
            if(exact_) {
                RVD_predicates::end_stats() ;
            }
        }

        unsigned int current_facet() const { return current_facet_ ; }
        Mesh* current_mesh() const { return current_mesh_ ; }

    protected:
        
        void clip_by_cell(
            const DelaunaySkeleton* skel, 
            unsigned int v, 
            Mesh*& ping, Mesh*& pong
        ) {
            for(unsigned int i=skel->star_begin(v); i<skel->star_end(v); i++) {
                unsigned int w = skel->neighbor(i) ;
                const plane3& P = skel->bisector(i) ;
                pong->clear() ;
                clip_by_plane(ping, pong, P, v, w) ;
                Mesh* tmp = ping; ping = pong; pong = tmp ;
            }
        }

        Sign side(
            const plane3& P, 
            const VertexEdge& v,
            unsigned int cell_in, unsigned int cell_out
        ) const {
            if(exact_) {
                vec3 p1 = delaunay_->vertex(cell_in) ;
                vec3 p2 = delaunay_->vertex(cell_out) ;
                return RVD_predicates::side(p1, p2, v, this) ; 
            } 
            return (P.side(v) > 0.0 ? POSITIVE : NEGATIVE) ; 
        }

        void clip_by_plane(
            Mesh* ping, Mesh* pong, 
            const plane3& P, 
            unsigned int cell_in, unsigned int cell_out
        ) {
            for(unsigned int f=0; f<ping->nb_facets(); f++) {

                bool crosses = false ;

                unsigned int last_i = ping->facet_end(f) - 1 ;
                const VertexEdge* last_v = &(ping->vertex(last_i)) ;
                Sign last_status = side(P, *last_v, cell_in, cell_out) ;

                for(
                    unsigned int i=ping->facet_begin(f); 
                    i<ping->facet_end(f); i++
                ) {
                    const VertexEdge* v = &(ping->vertex(i)) ;
                    Sign status = side(P, *v, cell_in, cell_out) ; 
                    
                    if(status == 0) { // Can only occur in exact mode
                        crosses = true ;
                        if(!pong->in_facet()) { pong->begin_facet() ; }
                        pong->add_vertex(*v) ;
                    } else {
                        if(status != last_status) {
                            VertexEdge I ;
                            if(symbolic_ || exact_) {
                                ping->intersect(I, *last_v, *v, P, cell_out) ;
                            } else {
                                ping->intersect(I, *last_v, *v, P) ;
                            }
                            if(!pong->in_facet()) { pong->begin_facet() ; }
                            pong->add_vertex(I) ;
                            if(status > 0) {
                                pong->top_vertex().set_edge(*last_v) ;
                            } else {
                                pong->top_vertex().set_flag(
                                    VertexEdge::INTERSECT
                                ) ;
                            }
                        }
                        if(status > 0) {
                            if(!pong->in_facet()) { pong->begin_facet() ; }
                            pong->add_vertex(*v) ;
                        } else {
                            crosses = true ;
                        }
                    }
                    last_v = v ; last_status = status ; last_i = i ; 
                }
                if(pong->in_facet())  { 
                    pong->end_facet() ;  
                    pong->facet_info(pong->top()).id = ping->facet_info(f).id ;
                } 
                if(crosses) {
                    push(ping->facet_info(f).id, cell_out) ;
                }
            }
        }

        void visit(unsigned int f, unsigned int cell) {
            visited_[f].push_back(cell) ;
            m->mark_facet(f) ;
        }

        bool is_visited(unsigned int f, unsigned int cell) {
            for(unsigned int i=0; i<visited_[f].size(); i++) {
                if(visited_[f][i] == cell) { return true ; }
            }
            return false ;
        }

        void push(
            unsigned int f, 
            unsigned int cell, bool check = true
        ) {
            if(check && is_visited(f, cell)) { return ; }
            visit(f, cell) ;
            S.push(RVDStackItem(f, cell)) ;
        }

    private:
        RVDStack S ;
        Mesh M1, M2 ;
        TopoCellSet* visited_ ;
        unsigned int visited_size_ ;
        Mesh* m ;
        Delaunay* delaunay_ ;
        bool symbolic_ ;
        bool exact_ ;
        unsigned int current_facet_ ;
        Mesh* current_mesh_ ;

    private:
        RestrictedVoronoiDiagram(const thisclass& rhs) ;
        thisclass& operator=(const thisclass& rhs) ;        
    } ;
}

#endif
