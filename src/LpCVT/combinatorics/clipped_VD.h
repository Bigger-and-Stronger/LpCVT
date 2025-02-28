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

#ifndef __GEEX_CVT_CLIPPED_VD__
#define __GEEX_CVT_CLIPPED_VD__

#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/symbolic_vertex.h>
#include <map>
#include <vector>

namespace Geex {

    /**
     * Given a Mesh and a Delaunay,
     * ClippedVoronoiDiagram computes the
     * intersections between the interior of
     * the surface and the cells of the Voronoi
     * diagram.
     *
     */
    class ClippedVoronoiDiagram {
    public:

        enum FaceType {CVD_NONE, CVD_RVD, CVD_WALL, CVD_CLOSE, CVD_INSIDE} ;
        FaceType current_face_type() const { return current_face_type_ ; }

        struct Vertex : public vec3Sym {
        public:
            Vertex() : visited(false) {}
            Vertex(
                const VertexEdge& v
            ) : vec3Sym(v), visited(false) {
            }
            Vertex(const vec3& v, const SymbolicVertex& k_in) : 
                vec3Sym(v,k_in), visited(false) { 
            }
            bool visited ;
        } ;

    protected:
        // Used to store the borders of the 
        // Restricted Voronoi Cells.
        typedef std::map<SymbolicVertex, Vertex> VertexMap ;

        // Used as an argument of RestrictedVoronoiDiagram::for_each_facet(),
        // stores the borders of the Restricted Voronoi Cells in an array
        // and applies the user-specified action to the facets of the restricted
        // voronoi cells.
        template <class ACTION> class GetCellsBoundary {
        public:
            GetCellsBoundary(
                std::vector<VertexMap>& cells, const ACTION& action
            ) : cells_(cells), action_(action) {}
            void operator()(unsigned int v, Mesh* M) const {
                if(action_.get_restricted_cells()) {
                    for(unsigned int f=0; f<M->nb_facets(); f++) {
                        for(unsigned int i1=M->facet_begin(f); 
                            i1<M->facet_end(f); i1++) {
                            if(M->vertex(i1).flags == VertexEdge::INTERSECT) {
                                unsigned int i2=M->next_around_facet(f,i1) ;
                                cells_[v][M->vertex(i1).sym] = Vertex(M->vertex(i2)) ;
                            }
                        }
                    }
                }
                if(action_.process_cell(v)) {
                    for(unsigned int f=0; f<M->nb_facets(); f++) {
                        action_.begin_facet(v, -1) ;
                        for(unsigned int i=M->facet_begin(f); i<M->facet_end(f); i++) {
                            action_.vertex(M->vertex(i)) ;
                        }
                        action_.end_facet() ;
                    }
                }
            }
        private:
            std::vector<VertexMap>& cells_ ;
            ACTION action_ ;
        } ;

        template <class TRIACTION> class TriangleAction {
        public:
            TriangleAction(
                const TRIACTION& action, bool recompute_restricted_cells = true
            ) : action_(action), recompute_rvc_(recompute_restricted_cells) {
            }
            bool get_restricted_cells() const             { return recompute_rvc_ ; }
            bool process_cell(unsigned int i) const       { return true ; }
            void begin_facet(int i, unsigned int j) const { vertices_.resize(0); i_ = i ; j_ = j ; }
            void vertex(const vec3Sym& v) const           { vertices_.push_back(v) ; }
            void end_facet() const {  
                if(vertices_.size() > 2) {
                    //   Facets on the boundary, i.e. from restricted Voronoi cells,
                    // are traversed in reverse orientation, therefore we invert the
                    // triangles when traversing them.
                    if(j_ == -1) {
                        for(unsigned int i=1; i+1<vertices_.size(); i++) {
                            action_(i_, j_, vertices_[0], vertices_[i], vertices_[i+1]) ;
                        }
                    } else {
                        for(unsigned int i=1; i+1<vertices_.size(); i++) {
                            action_(i_, j_, vertices_[0], vertices_[i+1], vertices_[i]) ;
                        }
                    }
                }
            }
        protected:
            const TRIACTION& action_ ;
            bool recompute_rvc_ ;
            mutable std::vector<vec3Sym> vertices_ ;
            mutable unsigned int i_ ;
            mutable int j_ ;
        } ;
        

    public:
        ClippedVoronoiDiagram(Delaunay* delaunay, Mesh* M) : RVD_(delaunay,M) {
            RVD_.set_symbolic(true) ;
            traverse_inner_cells_ = true ;
            geometry_ = true ;
        }

        RestrictedVoronoiDiagram& RVD() { return RVD_ ; }
        const RestrictedVoronoiDiagram& RVD() const { return RVD_ ; }
        const Delaunay* delaunay() const { return RVD_.delaunay() ; }
        Delaunay* delaunay() { return RVD_.delaunay() ; }

        void set_traverse_inner_cells(bool x) {
            traverse_inner_cells_ = x ;
        }

        void set_geometry(bool x) { 
            geometry_ = x ; 
        }

        /**
         * ACTION needs to define the following functions:
         * void operator()(
         *    unsigned int i, int j, const vec3Sym& v1, const vec3Sym& v2, const vec3Sym& v3
         * )
         */
        template <class TRIACTION> void for_each_triangle(
            const TRIACTION& action, bool recompute_rvc = true
        ) {
            for_each_facet(TriangleAction<TRIACTION>(action, recompute_rvc)) ;
        }

        /**
         * ACTION needs to define the following functions:
         * bool get_restricted_cells() -> returns true if restricted
         *                                cells need to be recomputed.
         * bool process_cell(unsigned int i) -> returns true if the cell
         *                                      associated with primal vertex i
         *                                      needs to be processed.
         * void begin_facet(unsigned int i, unsigned int j) -> starts the facet
         *                                      associated with bisector [i,j]
         *                                      j=-1 for facets on border.
         * void vertex(const vec3Sym& v)
         * void end_facet()
         */
        template <class ACTION> void for_each_facet(const ACTION& action) {
            current_face_type_ = CVD_NONE ;
            if(action.get_restricted_cells()) {
                restricted_cells_.resize(RVD_.delaunay()->nb_vertices()) ;
                for(unsigned int i=0; i<restricted_cells_.size(); i++) {
                    restricted_cells_[i].clear() ;
                }
                current_face_type_ = CVD_RVD ;
                RVD_.for_each_facet(GetCellsBoundary<ACTION>(restricted_cells_, action)) ;
            }
            if(traverse_inner_cells_) {
                cell_visited_.resize(RVD_.delaunay()->nb_vertices()) ;
                cell_on_stack_.resize(RVD_.delaunay()->nb_vertices()) ;
                std::fill(cell_visited_.begin(), cell_visited_.end(), false) ;
                std::fill(cell_on_stack_.begin(), cell_on_stack_.end(), false) ;
            }
            for(unsigned int iv=0; iv<RVD_.delaunay()->nb_vertices(); iv++) {
                traverse_cell_ = action.process_cell(iv) ;
                if(traverse_cell_ || traverse_inner_cells_) {
                    current_face_type_ = CVD_WALL ;
                    begin_cell(iv) ;
                    for(VertexMap::iterator it = RVC_->begin() ; it != RVC_->end(); it++) {
                        if(!it->second.visited) { 
                            get_contour(it, action) ; 
                        }
                    }       
                    for(unsigned int i=0; i<cell_.nb_facets(); i++) {
                        close_facet(i, action) ;
                        if(facet_vertices_[i].size() > 0) {
                            unsigned int jv = cell_.facet_bisector(i) ;
                            facet_visited_.insert(jv) ;
                            if(traverse_cell_) {
                                action.begin_facet(iv_, jv) ;
                                for(unsigned int k=0; k<facet_vertices_[i].size(); k++) {
                                    action.vertex(facet_vertices_[i][k]) ;
                                }
                                action.end_facet() ;
                            }
                        }
                    }
		    current_face_type_ = CVD_CLOSE ;
  		    propagate_to_adjacent_facets(action) ;
                    end_cell() ;
                    if(traverse_inner_cells_) {
                        cell_visited_[iv_] = (RVC_->size() > 0) ;
                    }
                }
            }
            if(traverse_inner_cells_) {
                current_face_type_ = CVD_INSIDE ;
                propagate_to_adjacent_cells(action) ;
            }
        }

    protected:

        template <class ACTION> void propagate_to_adjacent_cells(const ACTION& action) {
            while(!cell_S_.empty()) {
                iv_ = cell_S_.top() ;
                cell_S_.pop() ;
                if(
                    !cell_visited_[iv_] 
                ) {
                    cell_visited_[iv_] = true ;
                    if(action.process_cell(iv_)) {
                        RVD_.delaunay()->get_voronoi_cell(iv_, cell_, geometry_) ;
                        for(unsigned int fi=0; fi<cell_.nb_facets(); fi++) {
                            apply_action_to_cell_facet(fi, action) ;
                            if(cell_.facet_bisector(fi) >= 0) {
                                push_cell(cell_.facet_bisector(fi)) ;
                            }
                        }
                    }
                }
            }
        }

        template <class ACTION> void propagate_to_adjacent_facets(const ACTION& action) {
            while(!facet_S_.empty()) {
                unsigned int f = facet_S_.top() ; facet_S_.pop() ;
                if(facet_visited_.find(f) == facet_visited_.end()) {
                    facet_visited_.insert(f) ;
                    unsigned int fi = cell_.find_facet(f) ;
                    if(traverse_cell_) {
                        apply_action_to_cell_facet(fi, action) ;
                    }
                    if(traverse_inner_cells_) {
                        push_cell(f) ;
                    } 
                    for(unsigned int i = cell_.facet_begin(fi); i<cell_.facet_end(fi); i++) {
                        int g = cell_.edge_bisector(i) ;
                        if(g >= 0 && facet_visited_.find(g) == facet_visited_.end()) {
                            facet_S_.push(g) ;
                        }
                    }
                }
            }
        }

        template <class ACTION> void apply_action_to_cell_facet(unsigned int fi, const ACTION& action) {
            unsigned int f = cell_.facet_bisector(fi) ;
            action.begin_facet(iv_, f) ;
            for(unsigned int i = cell_.facet_begin(fi); i<cell_.facet_end(fi); i++) {
                SymbolicVertex k ; 
                k.insert(cell_.edge_bisector(i)+1) ;
                k.insert(cell_.edge_bisector(cell_.prev_around_facet(fi, i))+1) ;
                k.insert(f+1) ;
                if(geometry_) {
                    action.vertex(Vertex(cell_.vertex(i),k)) ;
                } else {
                    action.vertex(Vertex(vec3(0,0,0),k)) ;
                }
            }
            action.end_facet() ;
        }

        // initiates computation of cell associated with primal vertex iv.
        void begin_cell(unsigned int iv) {
            iv_ = iv ;
            RVC_ = &(restricted_cells_[iv]) ;
            RVD_.delaunay()->get_voronoi_cell(iv, cell_) ;
            facet_vertices_.resize(cell_.nb_facets()) ;
            for(unsigned int i=0; i<facet_vertices_.size(); i++) {
                facet_vertices_[i].resize(0) ;
            }
            facet_first_.resize(cell_.nb_facets()) ;
            facet_cur_.resize(cell_.nb_facets()) ;
            facet_enter_.resize(cell_.nb_facets()) ;
            facet_visited_.clear() ;
            std::fill(facet_first_.begin(), facet_first_.end(), -1) ;
            std::fill(facet_cur_.begin(), facet_cur_.end(), -1) ;
            std::fill(facet_enter_.begin(), facet_enter_.end(), -1) ;
        }

        void end_cell() {
            for(VertexMap::iterator it = RVC_->begin(); it != RVC_->end(); it++) {
                it->second.visited = false ;
            }
        }
        
        bool has_several_contours() {
            VertexMap::iterator it = RVC_->begin() ;
            unsigned int count = 0 ;
            do {
                count++ ;
                it = next(it) ;
            } while(it != RVC_->begin()) ;
            return (count < RVC_->size()) ;
        }

        unsigned int nb_contours() {
            unsigned int result = 0 ;
            for(VertexMap::iterator it = RVC_->begin(); it != RVC_->end(); it++) {
                if(!it->second.visited) {
                    result++ ;
                    VertexMap::iterator jt = it ;
                    do {
                        jt->second.visited = true ;
                        jt = next(jt) ;
                    } while(jt != it) ;
                }
            }
            for(VertexMap::iterator it = RVC_->begin(); it != RVC_->end(); it++) {
                it->second.visited = false ;
            }
            return result ;
        }


        // finds a contour (i.e. a closed loop) in the border of
        // the current Restricted Voronoi Cell.
        template <class ACTION> void get_contour(
            VertexMap::iterator from, const ACTION& action
        ) {
            from = contour_start(from) ;
            if(!is_corner(from)) {
                // Case 1: contour does not have corners, all in same bisector
                //   The contour is directly processed.
                unsigned int jv = from->second.sym[2]-1;
                facet_visited_.insert(jv) ;
                if(traverse_cell_) {
                    action.begin_facet(iv_, jv) ;
                    VertexMap::iterator it = from ;
                    do {
                        action.vertex(it->second) ;
                        it->second.visited = true ;
                        it = next(it) ;
                    } while(it != from) ;
                    action.end_facet() ;
                }
            } else {
                // Case 2: contour has corners
                //   The contour is stored.
                VertexMap::iterator it = from ;
                int cur_edge_bisector = edge_bisector(from) ;
                assert(cur_edge_bisector >= 0) ;
                unsigned int cur_facet = cell_.find_facet(cur_edge_bisector) ;
                enter_facet(cur_facet, it, action) ;
                do {
                    it->second.visited = true ;
                    // Detect whether we changed facet in the 
                    // current Voronoi cell.
                    unsigned int it_edge_bisector = edge_bisector(it) ;
                    if(it_edge_bisector != cur_edge_bisector) {
                        leave_facet(cur_facet, it, action) ;
                        cur_edge_bisector = it_edge_bisector ;
                        cur_facet = cell_.find_facet(cur_edge_bisector) ;
                        enter_facet(cur_facet, it, action) ;
                    } else {
                        facet_vertices_[cur_facet].push_back(it->second) ;
                    }
                    it = next(it) ;
                } while(it != from) ;
                leave_facet(cur_facet, it, action) ;
            }
        }

        template <class ACTION> 
        void enter_facet(unsigned int f, VertexMap::iterator it, const ACTION& action) {
            unsigned int eb = prev_edge_bisector(it) ;
            if(facet_first_[f] == -1) {
                facet_first_[f] = eb ;
                facet_enter_[f] = eb ;
            } else {
                connect(f, facet_cur_[f], eb, action) ;
                facet_enter_[f] = eb ;
            }
            facet_vertices_[f].push_back(it->second) ;
        }

        template <class ACTION>
        void leave_facet(unsigned int f, VertexMap::iterator it, const ACTION& action) {
            facet_cur_[f] = edge_bisector(it) ;
            facet_vertices_[f].push_back(it->second) ;
        }

        template <class ACTION>
        void close_facet(unsigned int f, const ACTION& action) {
            if(facet_first_[f] != -1) {
                assert(facet_cur_[f] != -1) ;
                connect(f, facet_cur_[f], facet_first_[f], action) ;
            } 
        }

        template <class ACTION>
        void facet_component(unsigned int f, const ACTION& action) {
            if(facet_vertices_[f].size() > 0) {
                unsigned int jv = cell_.facet_bisector(f) ;
                facet_visited_.insert(jv) ;
                if(traverse_cell_) {
                    action.begin_facet(iv_, jv) ;
                    for(unsigned int k=0; k<facet_vertices_[f].size(); k++) {
                        action.vertex(facet_vertices_[f][k]) ;
                    }
                    action.end_facet() ;
                }
            }
            facet_vertices_[f].resize(0) ;
            facet_first_[f] = -1 ;
            facet_cur_[f] = -1 ;
            facet_enter_[f] = -1 ;
        }

        // Insert the vertices of the Voronoi facet that are met between
        // bisectors eb1 and eb2 when turning around the facet.
        template <class ACTION>
        void connect(unsigned int f, unsigned int eb1, unsigned int eb2, const ACTION& action) {
            unsigned int eb0 = facet_enter_[f] ;
            if(eb1 != eb2 && eb1 != eb0) {
                unsigned int ivv = cell_.facet_begin(f) ;
                while(cell_.edge_bisector(ivv) != eb1) {
                    ivv++ ;
                    assert(ivv < cell_.facet_end(f)) ;
                }
                unsigned int f_bisector = cell_.facet_bisector(f) ;
                int eb ;
                do {
                    unsigned int prev_bisector = cell_.edge_bisector(ivv) ;
                    ivv = cell_.next_around_facet(f, ivv);
                    {
                        SymbolicVertex k ; 
                        k.insert(f_bisector+1) ;
                        k.insert(prev_bisector+1) ;
                        prev_bisector = cell_.edge_bisector(ivv) ;
                        k.insert(prev_bisector+1) ;
                        facet_vertices_[f].push_back(Vertex(cell_.vertex(ivv), k)) ; 
                    }
                    eb = cell_.edge_bisector(ivv) ;
                    if(eb != eb2 && eb != eb0) { 
                        if(eb != -1) {
                            facet_S_.push(eb) ; 
                        } 
                    }
                } while(eb != eb2 && eb != eb0) ;
                if(eb == eb0) {
                    facet_component(f, action) ;
                }
            } else if(eb1 == eb0) {
                facet_component(f, action) ;
            }
        }

        //   For contours with corners (i.e. vertices that are on two bisectors),
        // we systematically start from a corner (makes the code simpler !)
        //   Note that all contours do not necessarily have corners (e.g. tubular
        // zones that are undersampled).
        VertexMap::iterator contour_start(VertexMap::iterator from) {
            VertexMap::iterator it = from ;
            do {
                if(is_corner(it)) { return it ; }
                it = next(it) ;
            } while(it != from) ;
            return it ;
        }

        VertexMap::iterator next(VertexMap::iterator from) {
            VertexMap::iterator result = RVC_->find(from->second.sym) ;
            assert(result != RVC_->end()) ;        
            return result ;
        }

        // Corners are vertices that are on two bisectors.
        bool is_corner(VertexMap::iterator it) {
            const SymbolicVertex& s = it->second.sym ;
            return (s.size() == 3) && (s[1] > 0) && (s[2] > 0) ;
        }

        unsigned int edge_bisector(VertexMap::iterator it) {
            VertexMap::iterator jt = next(it) ;
            return edge_bisector(it->second.sym, jt->second.sym) ;
        }

        unsigned int prev_edge_bisector(VertexMap::iterator it) {
            return edge_bisector(it->first, it->second.sym) ;
        }

        // finds the bisector common to two vertices.
        static unsigned int edge_bisector(const SymbolicVertex& s1, const SymbolicVertex& s2) {
            SymbolicVertex I ;
            sets_intersect(s1, s2, I) ;
            // I is supposed to have one negative and one positive value,
            // the positive one is the last one
            assert(I.size() == 2) ;
            assert(I[0] < 0) ;
            assert(I[1] > 0) ;                            
            return I[1]-1 ;
        }

        void push_cell(unsigned int i) {
            if(!cell_on_stack_[i] && !cell_visited_[i]) {
                cell_on_stack_[i] = true ;
                cell_S_.push(i) ;
            }
        }

    private:
        unsigned int iv_ ;
        RestrictedVoronoiDiagram RVD_ ;
        bool traverse_cell_ ;
        VoronoiCell cell_ ;
        std::vector<VertexMap> restricted_cells_ ;
        VertexMap* RVC_ ;
        std::vector<int> facet_first_ ;
        std::vector<int> facet_cur_ ;
        std::vector<int> facet_enter_ ;
        std::vector< std::vector<Vertex> > facet_vertices_ ;
        std::set<unsigned int> facet_visited_ ;
        std::stack<unsigned int> facet_S_ ;
        bool traverse_inner_cells_ ;
        std::stack<unsigned int> cell_S_ ;
        std::vector<bool> cell_visited_ ;
        std::vector<bool> cell_on_stack_ ;
        bool geometry_ ;
        FaceType current_face_type_ ;
    } ;

}

#endif
