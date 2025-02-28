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

#ifndef __GEEX_CVT_SYMBOLIC_VERTEX__
#define __GEEX_CVT_SYMBOLIC_VERTEX__

#include <LpCVT/common/types.h>
#include <LpCVT/common/containers.h>

namespace Geex {

    typedef Numeric::uint64 SymbolicVertexKey ;

    typedef union {
        SymbolicVertexKey k ;
        Numeric::int32 words[2] ;
    } SymbolicVertexKeyAccess ;


    // This set of three integers encodes the equation of this vertex.
    // * Each positive entry i denotes the bisector of the segment that connects
    //   the center vertex to the i-th vertex (note that the center vertex
    //   needs to be stored elsewhere, but is known when a RVD or ClippedVD is used,
    //   since we know which dual cell we are processing).
    // * Each negative entry i denotes the i-th face in the boundary TriMesh.
    //   Note: indexing starts with 1 (resp. -1), 0 is kept for error codes.

    template <int N> class GenericSymbolicVertex : public small_set<int,N> {
    public:
        typedef small_set<int,N> baseclass ;
        typedef baseclass Sym ;
        // Returns a key formed of two integers packed into a single 64 bits integers.
        // They are used for chaining the boundary edges in a table (note that only
        // two equations are needed: since we just computed the intersection with a 
        // plane, all the boundary edges have this plane in their equation).
        SymbolicVertexKey key(int E) const {
            SymbolicVertexKeyAccess result ;
            if((*this)[0] == E) {
                result.words[0] = (*this)[1] ;
                result.words[1] = (*this)[2] ;
            } else if((*this)[1] == E) {
                result.words[0] = (*this)[0] ;
                result.words[1] = (*this)[2] ;
            } else {
                result.words[0] = (*this)[0] ;
                result.words[1] = (*this)[1] ;
            }
            return result.k ;
        }

        void add_bisector(int i) {
            baseclass::insert(i+1) ;
        }

        void add_boundary_facet(int i) {
            baseclass::insert(-i-1) ;
        }

        unsigned int nb_boundary_facets() const {
            unsigned int result = 0 ;
            for(
                typename baseclass::const_iterator it = baseclass::begin(); 
                it != baseclass::end() && *it<0; it++
            ) {
                result++ ;
            }
            return result ;
        }

        unsigned int nb_bisectors() const {
            unsigned int result = 0 ;
            for(
                typename baseclass::const_iterator it = baseclass::end() - 1; 
                it != baseclass::begin() - 1 && *it>0; it--
            ) {
                result++ ;
            }
            return result ;
        }

        unsigned int bisector(int i) const {
            return (baseclass::end()[-1-i])-1 ;
        }

        unsigned int boundary_facet(int i) const {
            return -(baseclass::begin()[i])-1 ;
        }

        bool has_bisector(int i) const {
            return baseclass::find(i+1) != baseclass::end() ;
        }

        bool has_boundary_facet(int i) const {
            return baseclass::find(-i-1) != baseclass::end() ;
        }
    } ;

    typedef GenericSymbolicVertex<3> SymbolicVertex ;

    inline bool operator<(const SymbolicVertex::Sym& k1, const SymbolicVertex::Sym& k2) { 
        if(k1[0] < k2[0]) { return true ;  }
        if(k1[0] > k2[0]) { return false ; }        
        if(k1[1] < k2[1]) { return true ;  }
        if(k1[1] > k2[1]) { return false ; }        
        if(k1[2] < k2[2]) { return true ;  }
        return false ;
    }

    inline bool operator==(
        const SymbolicVertex::Sym& k1, const SymbolicVertex::Sym& k2
    ) { 
        return 
            k1[0] == k2[0] &&
            k1[1] == k2[1] &&
            k1[2] == k2[2] ;
    }

    //__________________________________________________________________

    class vec3Sym : public vec3 {
    public:
        vec3Sym(
            const vec3& v, const SymbolicVertex& k
        ) : vec3(v), sym(k) {
        }
        vec3Sym(const vec3& v) : vec3(v) { }
        vec3Sym() { }
        SymbolicVertex sym ;
    } ;

}

#endif
