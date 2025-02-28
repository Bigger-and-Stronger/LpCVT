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

#ifndef __GEEX_CVT_DELAUNAY_SKEL__
#define __GEEX_CVT_DELAUNAY_SKEL__

#include <LpCVT/common/containers.h>
#include <LpCVT/common/types.h>
#include <stack>

namespace Geex {

    /**
     * DelaunaySkeleton stores the vertex graph of a Delaunay triangulation
     * and the equations of the bisectors associated with each primal edge.
     * All information is stored in a Compressed Row Storage array.
     */
    class DelaunaySkeleton {
    public:
        DelaunaySkeleton() {
            star_ptr_.push_back(0) ;
        }

        void clear() {
            // Note: we use resize(0) instead of clear()
            // since this is guaranteed to keep the reserved
            // memory (this is what we want since we keep 
            // clearing and filling the same DelaunaySkeleton).
            star_ptr_.resize(0) ; 
            neighbor_.resize(0) ; 
            bisector_.resize(0) ;
            star_ptr_.push_back(0) ;
        }

        unsigned int nb_vertices() const { return (unsigned int)star_ptr_.size() - 1 ; }
        unsigned int star_begin(unsigned int i) const {
            return star_ptr_[i] ;
        }
        unsigned int star_end(unsigned int i) const {
            return star_ptr_[i+1] ;
        }
        unsigned int nb_neighbors(unsigned int i) const {
            return star_end(i) - star_begin(i) ;
        }
        unsigned int neighbor(unsigned int j) const {
            return neighbor_[j] ;
        }
        const plane3& bisector(unsigned int j) const {
            return bisector_[j] ;
        }
        unsigned int neighbor(unsigned int i, unsigned int j) const {
            return neighbor(star_begin(i) + j) ;
        }
        const plane3& bisector(unsigned int i, unsigned int j) const {
            return bisector(star_begin(i) + j) ;
        }
        
        void begin_star() { }
        void add_to_star(unsigned int neighbor, const plane3& bisector) {
            neighbor_.push_back(neighbor) ;
            bisector_.push_back(bisector) ;
        }
        void end_star() {
            star_ptr_.push_back((unsigned int)neighbor_.size()) ;
        }

        // Used internally by copy_Delaunay_to_skel
        bool last_star_has_vertex(unsigned int i) const {
            for(
                unsigned int k = star_ptr_[star_ptr_.size() - 1]; k < neighbor_.size(); k++) {
                if(neighbor_[k] == i) { return true ; }
            }
            return false ;
        }

    private:
        std::vector<unsigned int> star_ptr_ ;
        std::vector<unsigned int> neighbor_ ;
        std::vector<plane3 > bisector_ ;
    } ;

//===========================================================================================


} 

#endif
