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


#include <LpCVT/algebra/F_Lp.h>
#include <LpCVT/algebra/voro_func.h>

#ifdef CVT_MULTITHREAD
#include <LpCVT/common/processor.h>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#endif

namespace Geex {
 
   
   //==========================================================================

    // Internal version: uses pointers instead of std::vectors for 'sym' and 'C'
    // (used by both monothread and multithread implementations)
    
    double compute_F_Lp_internal(
        bool volumic,                  // IN: false for surface meshing, true for volume meshing
        unsigned int p,                // IN: Lp norm to be used
        Mesh* mesh,                    // IN: the PLC
        unsigned int nb,               // IN: number of integration simplices
        const int* sym,                      // IN: 10 integers per integration simplex:
                                       //   - index of x0
                                       //   - symbolic representation of C1,C2,C3 (3 integers each)
        const vec3* C,                       // IN: C vertices of integration simplices (3 per integration simplex) 
        const std::vector<vec3>& X,    // IN: vertices
        const std::vector<plane3>& Q,  // IN: boundary facets
        const std::vector<mat3>& M,    // IN: anisotropy matrix (one per integration simplex)
        std::vector<double>& g         // OUT: gradient
    ) {
        assert(p >= 2 && p <= 16) ;
        assert((p/2)*2 == p) ;
        double f = 0.0 ;
        if(volumic) {
            switch(p) {
            case 2: {
                VoroFunc< IntegrationSimplex<2,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 4: {
                VoroFunc< IntegrationSimplex<4,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 6: {
                VoroFunc< IntegrationSimplex<6,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 8: {
                VoroFunc< IntegrationSimplex<8,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 10: {
                VoroFunc< IntegrationSimplex<10,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 12: {
                VoroFunc< IntegrationSimplex<12,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 14: {
                VoroFunc< IntegrationSimplex<14,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 16: {
                VoroFunc< IntegrationSimplex<16,TetVolume> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            }
        } else {
            switch(p) {
            case 2: {
                VoroFunc< IntegrationSimplex<2,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 4: {
                VoroFunc< IntegrationSimplex<4,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 6: {
                VoroFunc< IntegrationSimplex<6,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 8: {
                VoroFunc< IntegrationSimplex<8,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 10: {
                VoroFunc< IntegrationSimplex<10,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 12: {
                VoroFunc< IntegrationSimplex<12,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 14: {
                VoroFunc< IntegrationSimplex<14,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            case 16: {
                VoroFunc< IntegrationSimplex<16,TriArea> > F ;
                f = F.eval(mesh,nb,sym,C,X,Q,M,g) ;
            } break ;
            }
        }
        return f ;
    }
   

//=====================================================================


#ifdef CVT_MULTITHREAD

    /**
     * Stores all the parameters required to evaluate F_Lp on a set of integration
     * simplices. Manages also a local storage for the gradient.
     */
    class F_Lp_thread {
    public:
        void set_parameters(
            bool volumic_in,                  // IN: false for surface meshing, true for volume meshing
            unsigned int p_in,                // IN: Lp norm to be used
            Mesh* mesh_in,                    // IN: the PLC
            unsigned int nb_in,               // IN: number of integration simplices
            const int* sym_in,                // IN: 10 integers per integration simplex:
                                              //   - index of x0
                                              //   - symbolic representation of C1,C2,C3 (3 integers each)
            const vec3* C_in,                 // IN: C vertices of integration simplices (3 per integration simplex) 
            const std::vector<vec3>& X_in,    // IN: vertices
            const std::vector<plane3>& Q_in,  // IN: boundary facets
            const std::vector<mat3>& M_in     // IN: anisotropy matrix (one per integration simplex)
        ) {
            volumic = volumic_in ;
            p = p_in ;
            mesh = mesh_in ;
            nb = nb_in ;
            sym = sym_in ;
            C = C_in ;
            X = &X_in ;
            Q = &Q_in ;
            M = &M_in ;
            g.resize(X->size()*3) ;
        }

        void run() {
            std::fill(g.begin(), g.end(), 0.0) ;
            f = compute_F_Lp_internal(volumic, p, mesh, nb, sym, C, *X, *Q, *M, g) ;
        }

        bool volumic ;
        unsigned int p ;
        Mesh* mesh ;
        unsigned int nb ;
        const int* sym ;
        const vec3* C ;
        const std::vector<vec3>* X ;
        const std::vector<plane3>* Q ;
        const std::vector<mat3>* M ;
        double f ;
        std::vector<double> g ;
    } ;
    
    /**
     * User API function, 
     * multithread implementation.
     * 1) Partitions the set of integration simplices, 
     * 2) Calls compute_F_Lp_internal() in parallel,
     * 3) Gathers the result.
     */
    double compute_F_Lp(
        bool volumic,                  // IN: false for surface meshing, true for volume meshing
        unsigned int p,                // IN: Lp norm to be used
        Mesh* mesh,                    // IN: the PLC
        const std::vector<int>& sym,   // IN: 10 integers per integration simplex:
                                       //   - index of x0
                                       //   - symbolic representation of C1,C2,C3 (3 integers each)
        const std::vector<vec3>& C,    // IN: C vertices of integration simplices (3 per integration simplex) 
        const std::vector<vec3>& X,    // IN: vertices
        const std::vector<plane3>& Q,  // IN: boundary facets
        const std::vector<mat3>& M,    // IN: anisotropy matrix (one per integration simplex)
        std::vector<double>& g         // OUT: gradient
    ) {
        // Get number of cores from Processor class.
        unsigned int nb_threads = Processor::number_of_cores() ;
        std::cerr << "Evaluating F-Lp, using " << nb_threads << " threads " << std::endl ;
        std::vector<F_Lp_thread> threads(nb_threads) ;

        unsigned int nb_integration_simplices = (unsigned int)sym.size() / 10 ;
        unsigned int remaining = nb_integration_simplices ;
        unsigned int batch_size = nb_integration_simplices / nb_threads ;

        // Partition work (i.e. I and C arrays) into nb_threads blocs
        const int*  cur_I = &sym[0] ;
        const vec3* cur_C = &C[0] ;
        for(unsigned int i=0; i<nb_threads-1; i++) {
            threads[i].set_parameters(volumic, p, mesh, batch_size, cur_I, cur_C, X, Q, M) ;
            cur_I += batch_size * 10 ;
            cur_C += batch_size * 3 ;
            remaining -= batch_size ;
        }
        threads[threads.size()-1].set_parameters(volumic, p, mesh, remaining, cur_I, cur_C, X, Q, M) ;

        // Run the threads in parallel, using boost threads
        //   Note: Intel TBB has lower thread creation overhead (has a thread pool),
        // but we use here boost since LpCVT depends on CGAL that also depends on boost
        // (and for large meshes, thread creation time is negligible).
        boost::thread_group threads_impl ;
        for(unsigned int i=0; i<threads.size(); i++) {
            threads_impl.create_thread(boost::bind(&F_Lp_thread::run, &threads[i])) ;
        }
        threads_impl.join_all() ; // Waits for termination of all threads.

        // Gather the result
        double result = 0 ;
        std::fill(g.begin(), g.end(), 0.0) ;
        for(unsigned int i=0; i<threads.size(); i++) {
            result += threads[i].f ;
            for(unsigned int j=0; j<g.size(); j++) {
                g[j] += threads[i].g[j] ;
            }
        }
        return result ;
    }   

#else

    /**
     * User API function, 
     * monothread implementation.
     * This is just a wrapper around compute_F_Lp_internal()
     */
    double compute_F_Lp(
        bool volumic,                  // IN: false for surface meshing, true for volume meshing
        unsigned int p,                // IN: Lp norm to be used
        Mesh* mesh,                    // IN: the PLC
        const std::vector<int>& sym,   // IN: 10 integers per integration simplex:
                                       //   - index of x0
                                       //   - symbolic representation of C1,C2,C3 (3 integers each)
        const std::vector<vec3>& C,    // IN: C vertices of integration simplices (3 per integration simplex) 
        const std::vector<vec3>& X,    // IN: vertices
        const std::vector<plane3>& Q,  // IN: boundary facets
        const std::vector<mat3>& M,    // IN: anisotropy matrix (one per integration simplex)
        std::vector<double>& g         // OUT: gradient
    ) {
        unsigned int nb = sym.size() / 10 ;
        assert(sym.size() == nb*10) ;
        assert(C.size() == nb*3) ;
        return compute_F_Lp_internal(
            volumic, p, mesh, nb,
            &sym[0], &C[0], X, Q, M, g
        ) ;
    }   

#endif

}


