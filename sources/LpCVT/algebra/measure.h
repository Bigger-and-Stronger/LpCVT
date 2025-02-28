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

#ifndef __MEASURE_H__
#define __MEASURE_H__

#include <LpCVT/common/types.h>

//========================================================================================

namespace Geex {

    /**
     * Computes tetrahedron volume and its gradients.
     * See Appendix B.1 in the paper.
     * Note: the 1/6 factor is ignored.
     */
    class TetVolume {
    public:
        TetVolume() { }
        double operator()(
            const vec3& U1, const vec3& U2, const vec3& U3,
            vec3& dTdU1, vec3& dTdU2, vec3& dTdU3
        ) const {
            dTdU1 = cross(U2, U3) ;
            dTdU2 = cross(U3, U1) ;
            dTdU3 = cross(U1, U2) ;
            return dot(U1, dTdU1) ;
        }
    } ;

    /**
     * Computes triangle area and its gradients.
     * See Appendix B.1 in the paper.
     * Note: the 1/2 factor is ignored.
     */
    class TriArea {
    public:
        TriArea() { }
        double operator()(
            const vec3& U1, const vec3& U2, const vec3& U3,
            vec3& dTdU1, vec3& dTdU2, vec3& dTdU3
        ) const {
            vec3 N = cross(U1-U3,U2-U3) ;
            double T = length(N) ;
            if(::fabs(T) < 1e-10) {
                dTdU1 = vec3(0.0, 0.0, 0.0) ;
                dTdU2 = vec3(0.0, 0.0, 0.0) ;
                dTdU3 = vec3(0.0, 0.0, 0.0) ;
            } else {
                N = (1.0 / T)*N ;
                dTdU1 = cross(N,U3-U2);
                dTdU2 = cross(N,U1-U3);
                dTdU3 = cross(N,U2-U1);
            }
            return T ;
        }
    } ;

}

#endif
