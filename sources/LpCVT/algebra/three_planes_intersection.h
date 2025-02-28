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

#ifndef __THREE_PLANES_INTERSECTION_H__
#define __THREE_PLANES_INTERSECTION_H__

#include <LpCVT/common/types.h>

namespace Geex {

    /**
     * Computes Voronoi vertices and their derivatives.
     * (see Appendix B.2)
     */

    class ThreePlanesIntersection {
    public:
        ThreePlanesIntersection() { }

        // Configuration (D) (yellow vertices in Figure 4)
        void set(
            unsigned int i0_in, unsigned int i1_in, unsigned int i2_in, unsigned int i3_in
        ) {
            i0 = i0_in; i1 = i1_in ; i2 = i2_in ; i3 = i3_in ; nb_bisectors = 3 ;
            compute_W(X[i1]-X[i0], X[i2]-X[i0], X[i3]-X[i0]) ; 
            double p02 = X[i0].length2() ;
            compute_C(
                0.5*(X[i1].length2() - p02),
                0.5*(X[i2].length2() - p02),
                0.5*(X[i3].length2() - p02)
            ) ;
        }

        
        // Configuration (D) (yellow vertices in Figure 4), C given directly
        void set(
            unsigned int i0_in, unsigned int i1_in, unsigned int i2_in, unsigned int i3_in,
            const vec3& C_in
        ) {
            i0 = i0_in; i1 = i1_in ; i2 = i2_in ; i3 = i3_in ; nb_bisectors = 3 ;
            compute_W(X[i1]-X[i0], X[i2]-X[i0], X[i3]-X[i0]) ; 
            double p02 = X[i0].length2() ;
            C = C_in ;
        }

        
        // Configuration (C) (red vertices in Figure 4,6)
        void set(unsigned int i0_in, unsigned int i1_in, unsigned int i2_in, const plane3& Q1) {
            i0 = i0_in; i1 = i1_in ; i2 = i2_in ; nb_bisectors = 2 ;
            compute_W(X[i1]-X[i0], X[i2]-X[i0], Q1.normal()) ; 
            double p02 = X[i0].length2() ;
            compute_C(
                0.5*(X[i1].length2() - p02),
                0.5*(X[i2].length2() - p02),
                -Q1.d
            ) ;
        }

        // Configuration (C) (red vertices in Figure 4,6), C given directly
        void set(
            unsigned int i0_in, unsigned int i1_in, unsigned int i2_in, const plane3& Q1,
            const vec3& C_in
        ) {
            i0 = i0_in; i1 = i1_in ; i2 = i2_in ; nb_bisectors = 2 ;
            compute_W(X[i1]-X[i0], X[i2]-X[i0], Q1.normal()) ; 
            C = C_in ;
        }
        
        // Configuration (B) (green vertices in Figure 4,6)
        void set(unsigned int i0_in, unsigned int i1_in, const plane3& Q1, const plane3& Q2) {
            i0 = i0_in; i1 = i1_in ; nb_bisectors = 1 ;
            compute_W(X[i1]-X[i0], Q1.normal(), Q2.normal()) ; 
            double p02 = X[i0].length2() ;
            compute_C(
                0.5*(X[i1].length2() - p02),
                -Q1.d,
                -Q2.d
            ) ;
        }

        // Configuration (B) (green vertices in Figure 4,6), C given directly
        void set(
            unsigned int i0_in, unsigned int i1_in, const plane3& Q1, const plane3& Q2,
            const vec3& C_in
        ) {
            i0 = i0_in; i1 = i1_in ; nb_bisectors = 1 ;
            compute_W(X[i1]-X[i0], Q1.normal(), Q2.normal()) ; 
            C = C_in ;
        }


        // Configuration (A) (blue vertices in Figure 4,6)
        void set(const plane3& Q1, const plane3& Q2, const plane3& Q3) {
            nb_bisectors = 0 ;
            compute_W(Q1.normal(), Q2.normal(), Q3.normal()) ; 
            compute_C(
                -Q1.d,
                -Q2.d,
                -Q3.d
            ) ;
        }

        // Configuration (A) given directly (blue vertices in Figure 4,6)
        void set(const vec3& C_in) {
            nb_bisectors = 0 ;
            C = C_in ;
        }

        // Propagates the gradient from an integration simplex to the 
        // variables using the chain rule.
        void compose_gradient(const vec3& dFdC) {
            switch(nb_bisectors) {
            case 3: {
                compose_gradient_3_bisectors(dFdC) ;
            } break ;
            case 2: {
                compose_gradient_2_bisectors(dFdC) ;
            } break ;
            case 1: {
                compose_gradient_1_bisector(dFdC) ;
            } break ;
            }
        }

    protected:
        // The circumcenter C is given by MC = B
        // W0, W1, W2 are the columns of M^{-1}
        void compute_W(const vec3& N0, const vec3& N1, const vec3& N2) {
            W0 = cross(N1,N2) ;
            W1 = cross(N2,N0) ;
            W2 = cross(N0,N1) ;
            double Tinv = 1.0 / dot(N0,W0) ; 
            W0.x *= Tinv ;
            W0.y *= Tinv ;
            W0.z *= Tinv ;
            W1.x *= Tinv ;
            W1.y *= Tinv ;
            W1.z *= Tinv ;
            W2.x *= Tinv ;
            W2.y *= Tinv ;
            W2.z *= Tinv ;
        }

        void compute_C(double Dx, double Dy, double Dz) {
            C.x = Dx * W0.x + Dy * W1.x + Dz * W2.x ;
            C.y = Dx * W0.y + Dy * W1.y + Dz * W2.y ;
            C.z = Dx * W0.z + Dy * W1.z + Dz * W2.z ;
        }

        // gradient with configuration (D)
        // (see appendix B.2)
        void compose_gradient_3_bisectors(const vec3& dFdC) {

            mat3 J ; vec3 dFdCJ ;

            // ================  dC/dp0

            // dC/dx0
            double d = C.x - X[i0].x ;
            J(0,0) = d * (W0.x + W1.x + W2.x) ;
            J(1,0) = d * (W0.y + W1.y + W2.y) ;
            J(2,0) = d * (W0.z + W1.z + W2.z) ;

            // dC/dy0
            d = C.y - X[i0].y ;
            J(0,1) = d * (W0.x + W1.x + W2.x) ;
            J(1,1) = d * (W0.y + W1.y + W2.y) ;
            J(2,1) = d * (W0.z + W1.z + W2.z) ;
            
            // dC/dz0
            d = C.z - X[i0].z ;
            J(0,2) = d * (W0.x + W1.x + W2.x) ;
            J(1,2) = d * (W0.y + W1.y + W2.y) ;
            J(2,2) = d * (W0.z + W1.z + W2.z) ;
            
            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i0]   += dFdCJ.x ;
            g[3*i0+1] += dFdCJ.y ;
            g[3*i0+2] += dFdCJ.z ;

            //================= dC/dp1

            // dC/dx1
            d = X[i1].x - C.x ;
            J(0,0) = d * W0.x ;
            J(1,0) = d * W0.y ;        
            J(2,0) = d * W0.z ;        
            
            // dC/dy1
            d = X[i1].y - C.y ;
            J(0,1) = d * W0.x ;
            J(1,1) = d * W0.y ;        
            J(2,1) = d * W0.z ;        
            
            // dC/dz1
            d = X[i1].z - C.z ;
            J(0,2) = d * W0.x ;
            J(1,2) = d * W0.y ;        
            J(2,2) = d * W0.z ;        

            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i1]   += dFdCJ.x ;
            g[3*i1+1] += dFdCJ.y ;
            g[3*i1+2] += dFdCJ.z ;
            
            //================= dC/dp2

            // dC/dx2
            d = X[i2].x - C.x ;
            J(0,0) = d * W1.x ;
            J(1,0) = d * W1.y ;
            J(2,0) = d * W1.z ;

            // dC/dy2
            d = X[i2].y - C.y ;
            J(0,1) = d * W1.x ;
            J(1,1) = d * W1.y ;
            J(2,1) = d * W1.z ;

            // dC/dz2
            d = X[i2].z - C.z ;
            J(0,2) = d * W1.x ;
            J(1,2) = d * W1.y ;
            J(2,2) = d * W1.z ;

            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i2]   += dFdCJ.x ;
            g[3*i2+1] += dFdCJ.y ;
            g[3*i2+2] += dFdCJ.z ;
            

            //================= dC/dp3

            // dC/dx3
            d = X[i3].x - C.x ;
            J(0,0) = d * W2.x ;
            J(1,0) = d * W2.y ;
            J(2,0) = d * W2.z ;

            // dC/dy3
            d = X[i3].y - C.y ;
            J(0,1) = d * W2.x ;
            J(1,1) = d * W2.y ;
            J(2,1) = d * W2.z ;

            // dC/dz3
            d = X[i3].z - C.z ;
            J(0,2) = d * W2.x ;
            J(1,2) = d * W2.y ;
            J(2,2) = d * W2.z ;

            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i3]   += dFdCJ.x ;
            g[3*i3+1] += dFdCJ.y ;
            g[3*i3+2] += dFdCJ.z ;

        }


        // gradient with configuration (C)
        // (see appendix B.2)
        void compose_gradient_2_bisectors(const vec3& dFdC) {

            mat3 J ; vec3 dFdCJ ;

            // ================  dC/dp0

            // dC/dx0
            double d = C.x - X[i0].x ;
            J(0,0) = d * (W0.x + W1.x ) ;
            J(1,0) = d * (W0.y + W1.y ) ;
            J(2,0) = d * (W0.z + W1.z ) ;
            
            // dC/dy0
            d = C.y - X[i0].y ;
            J(0,1) = d * (W0.x + W1.x ) ;
            J(1,1) = d * (W0.y + W1.y ) ;
            J(2,1) = d * (W0.z + W1.z ) ;

            // dC/dz0
            d = C.z - X[i0].z ;
            J(0,2) = d * (W0.x + W1.x ) ;
            J(1,2) = d * (W0.y + W1.y ) ;
            J(2,2) = d * (W0.z + W1.z ) ;

            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i0]   += dFdCJ.x ;
            g[3*i0+1] += dFdCJ.y ;
            g[3*i0+2] += dFdCJ.z ;

            // ================  dC/dp1

            // dC/dx1
            d = X[i1].x - C.x ;
            J(0,0) = d * W0.x ;
            J(1,0) = d * W0.y ;        
            J(2,0) = d * W0.z ;        
            
            // dC/dy1
            d = X[i1].y - C.y ;
            J(0,1) = d * W0.x ;
            J(1,1) = d * W0.y ;        
            J(2,1) = d * W0.z ;        

            // dC/dz1
            d = X[i1].z - C.z ;
            J(0,2) = d * W0.x ;
            J(1,2) = d * W0.y ;        
            J(2,2) = d * W0.z ;        

            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i1]   += dFdCJ.x ;
            g[3*i1+1] += dFdCJ.y ;
            g[3*i1+2] += dFdCJ.z ;

            // ================  dC/dp2

            // dC/dx2
            d = X[i2].x - C.x ;
            J(0,0) = d * W1.x ;
            J(1,0) = d * W1.y ;
            J(2,0) = d * W1.z ;
            
            // dC/dy2
            d = X[i2].y - C.y ;
            J(0,1) = d * W1.x ;
            J(1,1) = d * W1.y ;
            J(2,1) = d * W1.z ;

            // dC/dz2
            d = X[i2].z - C.z ;
            J(0,2) = d * W1.x ;
            J(1,2) = d * W1.y ;
            J(2,2) = d * W1.z ;

            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i2]   += dFdCJ.x ;
            g[3*i2+1] += dFdCJ.y ;
            g[3*i2+2] += dFdCJ.z ;
        }

        // gradient with configuration (B)
        // (see appendix B.2)
        void compose_gradient_1_bisector(const vec3& dFdC) {

            mat3 J ; vec3 dFdCJ ;
            
            //============= dC/dp0
            
            // dC/dx0
            double d = C.x - X[i0].x ;
            J(0,0) = d * W0.x ;
            J(1,0) = d * W0.y ;
            J(2,0) = d * W0.z ;
            
            // dC/dy0
            d = C.y - X[i0].y ;
            J(0,1) = d * W0.x ;
            J(1,1) = d * W0.y ;
            J(2,1) = d * W0.z ;

            // dC/dz0
            d = C.z - X[i0].z ;
            J(0,2) = d * W0.x ;
            J(1,2) = d * W0.y ;
            J(2,2) = d * W0.z ;

            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i0]   += dFdCJ.x ;
            g[3*i0+1] += dFdCJ.y ;
            g[3*i0+2] += dFdCJ.z ;

            //============= dC/dp1

            // dC/dx1
            d = X[i1].x - C.x ;
            J(0,0) = d * W0.x ;
            J(1,0) = d * W0.y ;        
            J(2,0) = d * W0.z ;        

            // dC/dy1
            d = X[i1].y - C.y ;
            J(0,1) = d * W0.x ;
            J(1,1) = d * W0.y ;        
            J(2,1) = d * W0.z ;        
            
            // dC/dz1
            d = X[i1].z - C.z ;
            J(0,2) = d * W0.x ;
            J(1,2) = d * W0.y ;        
            J(2,2) = d * W0.z ;        
            
            matTvecmul(J, dFdC, dFdCJ) ;
            g[3*i1]   += dFdCJ.x ;
            g[3*i1+1] += dFdCJ.y ;
            g[3*i1+2] += dFdCJ.z ;
        }


    public:
        const vec3* X ;
        double* g ;
        unsigned int nb_bisectors ;
        unsigned int i0 ;
        unsigned int i1 ;
        unsigned int i2 ;
        unsigned int i3 ;
        vec3 W0,W1,W2 ;
        vec3 C ;
    } ;

}

#endif

