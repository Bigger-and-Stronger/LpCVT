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

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/delaunay_CGAL.h>


namespace Geex {

    Delaunay::Delaunay() : skeleton_dirty_(false), nb_vertices_(0) {
        vertices_ = 0 ;
    }
    
    Delaunay::~Delaunay() { }

    Delaunay* Delaunay::create(const std::string& name) {
		std::cerr << "Creating Delaunay implementation, using " << name << std::endl ;
        if(name == "CGAL") {
            return new Delaunay_CGAL ;
        } else {
            std::cerr << name << ": no such Delaunay implementation." << std::endl ;
            return 0 ;
        }
    }

    void Delaunay::get_tetras(std::vector<int>& tetras, bool finite_only) {
        // Not implemented in baseclass
        assert(0) ;
        return ;
    }

    void  Delaunay::get_facets(std::vector<int>& facets, bool finite_only) {
        // Not implemented in baseclass
        assert(0) ;
        return ;
    }
    
    void Delaunay::get_voronoi_cell(unsigned int i, VoronoiCell& cell, bool geometry) {
        // Not implemented in baseclass
        assert(0) ;
        return ;
    }


//=======================================================================

}
