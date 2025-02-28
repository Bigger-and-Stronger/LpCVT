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

#ifndef __GEEX_MATHEMATICS_TYPES__
#define __GEEX_MATHEMATICS_TYPES__

#include <LpCVT/common/vecg.h>
#include <LpCVT/common/matrix.h>
#include <LpCVT/common/plane.h>

#include <math.h>
#include <assert.h>

namespace Geex {

    typedef vec2g<double> vec2 ;
    typedef vec3g<double> vec3 ;
    typedef vec4g<double> vec4 ;
    
    typedef Matrix<double, 3> mat3 ;
    typedef Matrix<double, 4> mat4 ;

    typedef Plane<double> plane3 ;

    namespace Numeric {
#ifdef WIN32
        typedef unsigned __int64 uint64 ;
#else
        typedef unsigned long long int uint64 ;
#endif
        typedef int int32 ;
    }

    enum Sign { NEGATIVE=-1, ZERO=0, POSITIVE=1 } ;

    template <class T> inline Sign gx_sgn(T x) {
        return (x > 0) ? POSITIVE : (
            (x < 0) ? NEGATIVE : ZERO
        );
    }

    //================ Some utility functions ======================

    inline double mixed_product(const vec3& v1, const vec3& v2, const vec3& v3) {
        return dot(cross(v1,v2),v3) ;
    }

    // V <- M U
    inline void matvecmul(
        const mat3& M,
        const vec3& U,
        vec3& V
    ) {
        V.x = M(0,0) * U.x + M(0,1) * U.y + M(0,2) * U.z ;
        V.y = M(1,0) * U.x + M(1,1) * U.y + M(1,2) * U.z ;
        V.z = M(2,0) * U.x + M(2,1) * U.y + M(2,2) * U.z ;
    }

    // V <- M^T U
    inline void matTvecmul(
        const mat3& M,
        const vec3& U,
        vec3& V
    ) {
        V.x = M(0,0) * U.x + M(1,0) * U.y + M(2,0) * U.z ;
        V.y = M(0,1) * U.x + M(1,1) * U.y + M(2,1) * U.z ;
        V.z = M(0,2) * U.x + M(1,2) * U.y + M(2,2) * U.z ;
    }

    // to <- p1 * p2
    inline void vecmul(const vec3& p1, const vec3& p2, vec3& to) {
        to.x = p1.x*p2.x ;
        to.y = p1.y*p2.y ;
        to.z = p1.z*p2.z ;
    }

    // to <- p1 * p2 * p3
    inline void vecmul(const vec3& p1, const vec3& p2, const vec3& p3, vec3& to) {
        to.x = p1.x*p2.x*p3.x ;
        to.y = p1.y*p2.y*p3.y ;
        to.z = p1.z*p2.z*p3.z ;
    }

    // to += s.(p1*p2*p3)
    inline void vecmadd(double s, const vec3& p1, const vec3& p2, const vec3& p3, vec3& to) {
        to.x += s * p1.x * p2.x * p3.x ;
        to.y += s * p1.y * p2.y * p3.y ;
        to.z += s * p1.z * p2.z * p3.z ;
    }

    // to += s.p1 + t.p2
    inline void vecmadd(double s, const vec3& p1, double t, const vec3& p2, vec3& to) {
        to.x = s*p1.x + t*p2.x ;
        to.y = s*p1.y + t*p2.y ;
        to.z = s*p1.z + t*p2.z ;
    }

    inline double vecbar(const vec3& p) {
        return p.x + p.y + p.z ;
    }


}


#endif

