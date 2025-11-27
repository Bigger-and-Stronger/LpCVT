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

#ifndef __GEEX_CVT_VORONOI_CELL__
#define __GEEX_CVT_VORONOI_CELL__

#include <LpCVT/common/types.h>
#include <vector>

namespace Geex {

    /**
     * VoronoiCell stores the dual facets in a Compressed Row Storage array.
     * - Each facet knows the bisector it is on, and the list of vertices/edges.
     *    - Each vertex knows the tet it is dual to.
     *    - Each edge knows the other bisector it is on (an edge is defined as the
     * intersection between the facet bisector and the edge bisector).
     */
    class VoronoiCell {
    public:
        VoronoiCell() { facet_ptr_.push_back(0) ; }
        void clear() {
            facet_ptr_.resize(0) ;
            facet_bisector_.resize(0) ;
            edge_bisector_.resize(0) ;
            vertex_.resize(0) ;
            infinite_.resize(0) ;
            facet_ptr_.push_back(0) ;
        }

        unsigned int nb_facets() const { return (unsigned int)facet_ptr_.size() - 1 ; }

        unsigned int facet_begin(unsigned int f) const {
            return facet_ptr_[f] ;
        }

        unsigned int facet_end(unsigned int f) const {
            return facet_ptr_[f+1] ;
        }

        unsigned int nb_vertices(unsigned int f) const {
            return facet_end(f) - facet_begin(f) ;
        }

        unsigned int next_around_facet(unsigned int f, unsigned int i) const {
            return (i+1 == facet_end(f) ? facet_begin(f) : i+1) ;
        }

        unsigned int prev_around_facet(unsigned int f, unsigned int i) const {
            return (i == facet_begin(f) ? facet_end(f)-1 : i-1) ;
        }

        unsigned int facet_bisector(unsigned int f) const {
            return facet_bisector_[f] ;
        }

        int edge_bisector(unsigned int i) const {
            return edge_bisector_[i] ;
        }

        const vec3& vertex(unsigned int i) const {
            return vertex_[i] ;
        }

        bool vertex_is_infinite(unsigned int i) const {
            return infinite_[i] ;
        }

        void begin_facet(unsigned int f_bisector) {
            facet_bisector_.push_back(f_bisector) ;
        }

        void add_to_facet(
            int e_bisector, const vec3& v, bool infinite
        ) {
            edge_bisector_.push_back(e_bisector) ;
            vertex_.push_back(v) ;
            infinite_.push_back(infinite) ;
        }

        void add_to_facet(
            int e_bisector, bool infinite
        ) {
            edge_bisector_.push_back(e_bisector) ;
            infinite_.push_back(infinite) ;
        }

        void end_facet() {
            facet_ptr_.push_back((unsigned int)edge_bisector_.size()) ;            
        }

        unsigned int find_facet(unsigned int bisector) {
            for(unsigned int i=0; i<facet_bisector_.size(); i++) {
                if(facet_bisector_[i] == bisector) {
                    return i ;
                }
            }
            std::cerr << "bisector = " << bisector ;
            std::cerr << "facet = [" ;
            for(unsigned int i=0; i<facet_bisector_.size(); i++) {
                std::cerr << facet_bisector_[i] << " " ;
            }
            std::cerr << "]" << std::endl ;
            assert(0) ;
            return 0 ;
        }

    private:
        std::vector<unsigned int> facet_ptr_ ;
        std::vector<unsigned int> facet_bisector_ ;
        std::vector<int> edge_bisector_ ;
        std::vector<vec3> vertex_ ;
        std::vector<bool> infinite_ ; 
    } ;

}

#endif
