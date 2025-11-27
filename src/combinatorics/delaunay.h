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

#ifndef __GEEX_CVT_DELAUNAY__
#define __GEEX_CVT_DELAUNAY__

#include <LpCVT/combinatorics/delaunay_skel.h>
#include <LpCVT/combinatorics/voronoi_cell.h>

#include <string>

namespace Geex {

    /**
     * Abstract API for Delaunay triangulation in 3D.
     * Input: set of vertices. Output: Delaunay skeleton.
     */
    class Delaunay {
    public:
        Delaunay() ;
        virtual ~Delaunay() ;

        virtual void set_vertices(unsigned int n, const double* xyz) = 0 ;
        void set_vertices(unsigned int n, const vec3* V) { 
            set_vertices(n, &(V[0].x)) ;   
        }       
        void set_vertices(const std::vector<vec3>& V)    {
            set_vertices((unsigned int)V.size(), &(V[0])) ; 
        }

        unsigned int nb_vertices() const { return nb_vertices_ ; }
        const double* vertices() const { return vertices_ ; }
        vec3 vertex(unsigned int i) const { 
            return vec3(vertices_[3*i], vertices_[3*i+1], vertices_[3*i+2]); 
        }

        virtual unsigned int nearest_vertex_id(double x, double y, double z) const = 0 ;
        unsigned int nearest_vertex_id(const vec3& p) const { 
            return nearest_vertex_id(p.x, p.y, p.z) ; 
        }

        const DelaunaySkeleton* skeleton() const {
            if(skeleton_dirty_) { 
                update_skeleton() ; 
            }
            return &skeleton_ ;
        }

        /**
         * Entry point to factory.
         * For the moment, name = "CGAL" is supported.
         */
        static Delaunay* create(const std::string& name) ;

        virtual void get_tetras(std::vector<int>& tetras, bool finite_only) ;
        virtual void get_facets(std::vector<int>& facets, bool finite_only) ;
        virtual void get_voronoi_cell(unsigned int i, VoronoiCell& cell, bool geometry=true) ;

    protected:
        virtual void update_skeleton() const = 0 ;

    protected:
        mutable bool skeleton_dirty_ ;
        mutable DelaunaySkeleton skeleton_ ;
        const double* vertices_ ;
        unsigned int nb_vertices_ ;
    } ;


    //===========================================================================================

}

#endif

