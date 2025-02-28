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

#ifndef __VORO_FUNC_H__
#define __VORO_FUNC_H__

#include <LpCVT/algebra/three_planes_intersection.h>
#include <LpCVT/algebra/integration_simplex.h>
#include <LpCVT/combinatorics/mesh.h>

namespace Geex {

    /**
     * Computes F_{L_p} and its gradient from the symbolic information
     * and parameters.
     * Template parameter: 
     *    FUNC: IntegrationSimplex formula (see integration_simplex.h)
     * Note: VoroFunc keeps a reference to the original mesh, to solve
     * a problematic configuration: intersection between a bisector and
     * an edge shared by two coplanar facets on the boundary (occurs
     * frequently with CAD models).
     */

    template <class FUNC> class VoroFunc {
    public:
        VoroFunc() {
        }

        void set_intersection(
            int i0, int i1, int i2, int i3,
            ThreePlanesIntersection& I,
            const vec3& C
        ) {

            if(i1 >= 0) {
                if(i2 >= 0) {
                    if(i3 >= 0) {
                        I.set(i0, i1-1, i2-1, i3-1, C) ;
                    } else {
                        I.set(i0, i1-1, i2-1, Q[-i3-1], C) ;
                    }
                } else {
                    // Intersection between a bisector and an edge of the input mesh.
                    // To resist coplanar facets, the following code replaces the second
                    // facet with the plane orthogonal to the first one that passes through
                    // the common edge.
                    int q1 = -i2-1 ;
                    int q2 = -i3-1 ;
                    const plane3& Q1 = Q[q1] ;
                    vec3 v1,v2 ;
                    mesh->find_edge_extremities(q2, q1, v1, v2) ;
                    vec3 N = cross(v2 - v1, Q1.normal()) ;
                    plane3 Q2(N.x, N.y, N.z, -dot(N,v1)) ;
                    I.set(i0, i1-1, Q1, Q2, C) ;
                    // This replaces the following line, not resistant to coplanar facets:
                    // I.set(i0, i1-1, Q[-i2-1], Q[-i3-1], C) ;
                }
            } else {
                I.set(C) ;
            }
        }

        double eval(
            Mesh* mesh_in,                   // IN: the PLC 
            unsigned int nb,                 // IN: number of integration simplices
            const int*  I,                   // IN: 10 integers per integration simplex:
                                             //   - index of x0
                                             //   - symbolic representation of C1,C2,C3 (3 integers each)
            const vec3* C,                   // IN: geometric representation of C1,C2,C3 (3 per integration simplex)
            const std::vector<vec3>& X,      // IN: vertices
            const std::vector<plane3>& Qin,  // IN: boundary facets
            const std::vector<mat3>& M,      // IN: anisotropy matrix (one per integration simplex)
            std::vector<double>& gradient    // OUT: gradient
        ) {
            mesh = mesh_in ;
            ThreePlanesIntersection Ix[3] ;

            double result = 0.0 ;
            Q = &(Qin[0]) ;
            g = &(gradient[0]) ;

            for(unsigned int i=0; i<3; i++) {
                Ix[i].g = g ;
                Ix[i].X = &(X[0]) ;
            }

            vec3 dFdp0, dFdC[3] ;
            vec3 U,V,W ; mat3 Aniso ;

            int prev_i0 = -1 ;
            int prev_f  = -1 ;
            
            for(unsigned int i=0; i<nb; i++) {
                set_intersection(I[0], I[1], I[2], I[3], Ix[0], C[0]) ;
                set_intersection(I[0], I[4], I[5], I[6], Ix[1], C[1]) ;
                set_intersection(I[0], I[7], I[8], I[9], Ix[2], C[2]) ;

                result += F.eval(
                    X[I[0]], Ix[0].C, Ix[1].C, Ix[2].C,
                    M[i],
                    dFdp0, dFdC[0], dFdC[1], dFdC[2]
                ) ;

                Ix[0].compose_gradient(dFdC[0]) ;
                Ix[1].compose_gradient(dFdC[1]) ;
                Ix[2].compose_gradient(dFdC[2]) ;

                // Compose grad for p0
                g[3*I[0]]   += dFdp0.x ;
                g[3*I[0]+1] += dFdp0.y ;
                g[3*I[0]+2] += dFdp0.z ;

                I += 10 ;
                C += 3 ;
            }
            
            return result ;
        }
    private:
        const vec3* C ;
        const plane3* Q ;
        double* g ;
        FUNC F ;
        Mesh* mesh ;
    } ;

}

#endif
