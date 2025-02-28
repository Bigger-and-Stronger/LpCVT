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

#ifndef __INTEGRATION_SIMPLEX_H__
#define __INTEGRATION_SIMPLEX_H__

#include <LpCVT/algebra/measure.h>

namespace Geex {

    /**
     * Computes F_{L_p}^T and its gradients.
     * See Appendices A and B.1 in the paper.
     * Template parameters:
     *    P       : integer (for specifying the norm L_P)
     *    MEASURE : TriArea or TetVolume (see measure.h)
     * Rem: the binomial coefficient (n+p \choose n) is ignored
     */

    template <int P, class MEASURE> class IntegrationSimplex {
    public:
        enum { degree = P, nb_coeffs = ((P+1)*(P+2))/2, nb_dcoeffs = nb_coeffs - (P+1) } ;

        IntegrationSimplex() {
            // Pre-compute indices for evaluating F_{L_p}^T
            // (see Appendix A)
            {
                unsigned int cur=0;
                for(unsigned int alpha=0; alpha<=P; alpha++) {
                    for(unsigned int beta=0; beta<=P-alpha; beta++) {
                        unsigned int gamma = P-alpha-beta ;
                        E_pow[cur][0] = alpha ;
                        E_pow[cur][1] = beta ;
                        E_pow[cur][2] = gamma ;
                        cur++ ;
                    }
                }
            }
            // Pre-compute indices for evaluating \nabla F_{L_p}^T
            // (see Appendix B.1)
            {
                unsigned int cur_dU1=0, cur_dU2=0, cur_dU3=0 ;
                for(unsigned int alpha=0; alpha<=P; alpha++) {
                    for(unsigned int beta=0; beta<=P-alpha; beta++) {
                        unsigned int gamma = P-alpha-beta ;
                        if(alpha != 0) {
                            dE_pow[0][cur_dU1][0] = alpha ;
                            dE_pow[0][cur_dU1][1] = beta ;
                            dE_pow[0][cur_dU1][2] = gamma ;
                            cur_dU1++ ;
                        }
                        if(beta != 0) {
                            dE_pow[1][cur_dU2][0] = alpha ;
                            dE_pow[1][cur_dU2][1] = beta ;
                            dE_pow[1][cur_dU2][2] = gamma ;
                            cur_dU2++ ;
                        }
                        if(gamma != 0) {
                            dE_pow[2][cur_dU3][0] = alpha ;
                            dE_pow[2][cur_dU3][1] = beta ;
                            dE_pow[2][cur_dU3][2] = gamma ;
                            cur_dU3++ ;
                        }
                    }
                }                
            }
            U_pow[0][0] = vec3(1.0, 1.0, 1.0) ;
            U_pow[1][0] = vec3(1.0, 1.0, 1.0) ;
            U_pow[2][0] = vec3(1.0, 1.0, 1.0) ;
        }

        double eval(
            const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3,
            const mat3& M,
            vec3& dFdp0, vec3& dFdp1, vec3& dFdp2, vec3& dFdp3
        ) {
            {
                vec3 U = p1 - p0 ;
                matvecmul(M, U, U_pow[0][1]) ;
                U=p2-p0 ;
                matvecmul(M, U, U_pow[1][1]) ;
                U=p3-p0 ;
                matvecmul(M, U, U_pow[2][1]) ;
            }

            for(unsigned int i=2; i<=P; i++) {
                vecmul(U_pow[0][1], U_pow[0][i-1], U_pow[0][i]) ;
                vecmul(U_pow[1][1], U_pow[1][i-1], U_pow[1][i]) ;
                vecmul(U_pow[2][1], U_pow[2][i-1], U_pow[2][i]) ;
            }

            // Computation of function value.
            double E = 0.0 ;
            for(unsigned int i=0; i<nb_coeffs; i++) {
                vec3 W ;
                unsigned int alpha = E_pow[i][0] ;
                unsigned int beta  = E_pow[i][1] ;
                unsigned int gamma = E_pow[i][2] ;
                vecmul(U_pow[0][alpha], U_pow[1][beta], U_pow[2][gamma], W) ;
                E += vecbar(W) ;
            }

            // Computation of gradient
            vec3 dEdU1(0,0,0), dEdU2(0,0,0), dEdU3(0,0,0) ;
            for(unsigned int i=0; i<nb_dcoeffs; i++) {
                {
                    unsigned int alpha = dE_pow[0][i][0] ;
                    unsigned int beta  = dE_pow[0][i][1] ;
                    unsigned int gamma = dE_pow[0][i][2] ;
                    vecmadd(alpha, U_pow[0][alpha-1], U_pow[1][beta], U_pow[2][gamma], dEdU1) ;
                }
                {
                    unsigned int alpha = dE_pow[1][i][0] ;
                    unsigned int beta  = dE_pow[1][i][1] ;
                    unsigned int gamma = dE_pow[1][i][2] ;
                    vecmadd(beta, U_pow[0][alpha], U_pow[1][beta-1], U_pow[2][gamma], dEdU2) ;
                }
                {
                    unsigned int alpha = dE_pow[2][i][0] ;
                    unsigned int beta  = dE_pow[2][i][1] ;
                    unsigned int gamma = dE_pow[2][i][2] ;
                    vecmadd(gamma, U_pow[0][alpha], U_pow[1][beta], U_pow[2][gamma-1], dEdU3) ;
                }
            }


            // Compute tet measure and its 
            // derivatives relative to U1, U2 and U3. 
            vec3 dTdU1, dTdU2, dTdU3 ;
            double T = measure_(
                U_pow[0][1], U_pow[1][1], U_pow[2][1],
                dTdU1, dTdU2, dTdU3
            ) ;

            // Assemble dF = E.d|T| + |T|.dE
            // Rem: anisotropy matrix needs to be transposed
            // grad(F(MX)) = J(F(MX))^T = (gradF^T M)^T = M^T grad F

            vec3 dFdU1, dFdU2, dFdU3 ;
            vecmadd(E, dTdU1, T, dEdU1, dFdU1) ;
            matTvecmul(M, dFdU1, dFdp1) ;
            vecmadd(E, dTdU2, T, dEdU2, dFdU2) ;
            matTvecmul(M, dFdU2, dFdp2) ;
            vecmadd(E, dTdU3, T, dEdU3, dFdU3) ;
            matTvecmul(M, dFdU3, dFdp3) ;

            // Gradient relative to p1 (resp. p2,p3) = gradient relative to U1 (resp U2,U3)
            // Gradient relative to p0 is equal to minus gradient relative to U1,U2 and U3
            dFdp0.x = -dFdp1.x-dFdp2.x-dFdp3.x ;
            dFdp0.y = -dFdp1.y-dFdp2.y-dFdp3.y ;
            dFdp0.z = -dFdp1.z-dFdp2.z-dFdp3.z ;

            return T*E ;
        }

    private:
        unsigned int E_pow[nb_coeffs][3] ;
        unsigned int dE_pow[3][nb_dcoeffs][3] ;
        vec3 U_pow[3][P+1] ;
        MEASURE measure_ ;
    } ;

}

#endif
