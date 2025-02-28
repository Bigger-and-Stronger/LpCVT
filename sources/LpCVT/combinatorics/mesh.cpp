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


#include <LpCVT/combinatorics/mesh.h>
#include <LpCVT/common/line_stream.h>
#include <map>
#include <fstream>

#ifdef WITH_METIS
extern "C" {
#include <METIS/metis.h>
}
#endif

namespace Geex {

    static bool has_edge(
        const std::vector<unsigned int>& facet_vertex_index, 
        unsigned int facet_begin, unsigned int facet_end,
        unsigned int v1, unsigned int v2
    ) {
        unsigned int facet_n = facet_end - facet_begin ;
        for(unsigned int i=0; i<facet_n; i++) {
            unsigned int w1 = facet_vertex_index[facet_begin + i] ;
            unsigned int w2 = facet_vertex_index[facet_begin + ((i + 1) % facet_n)] ;
            if(w1 == v1 && w2 == v2) {
                return true ;
            }
        }
        return false ;
    }


    unsigned int Mesh::load(const std::string& filename) {
        std::cerr << "Mesh: Loading " << filename << std::endl ;
        std::ifstream input(filename.c_str()) ;
        if(!input) {
            std::cerr << "could not open file." << std::endl ;
            return 0 ;
        }
        clear() ;
        LineInputStream in(input) ;
        std::vector<vec3> vertex ;
        std::vector< std::vector<unsigned int> > star ;

        // Step 1: load vertices and triangles
        // (and keep track of stars, i.e. lists
        //  of triangles incident to each vertex).
        unsigned int v_offset = 0 ;
        while(!in.eof()) {
            in.get_line() ;
            std::string keyword ;
            
            in >> keyword ;
            
            if(keyword == "v") {
                vec3 p ;
                in >> p ;
                vertex.push_back(p) ;
                star.push_back( std::vector<unsigned int>() ) ;
            } else if(keyword == "HEADER") {
                v_offset = (unsigned int)vertex.size() ;
            } else if(keyword == "VRTX" || keyword == "PVRTX" || keyword == "ATOM") {
                int idx ;
                vec3 p ;
                in >> idx >> p ;
                vertex.push_back(p) ;
                star.push_back( std::vector<unsigned int>() ) ;
            } else if(keyword == "TRGL") {
                int i1,i2,i3 ;
                in >> i1 >> i2 >> i3 ;
                std::vector<unsigned int> cur_facet ;
                cur_facet.push_back(i1-1+v_offset) ;
                cur_facet.push_back(i2-1+v_offset) ;
                cur_facet.push_back(i3-1+v_offset) ;
                unsigned int f = nb_facets() ;
                begin_facet() ;
                for(unsigned int i=0; i<cur_facet.size(); i++) {
                    unsigned int v = cur_facet[i] ;
                    add_vertex(VertexEdge(vertex[v])) ;
                    top_vertex().set_flag(VertexEdge::ORIGINAL) ;
                    vertex_index_.push_back(v) ;
                    star[v].push_back(f) ;
                }
                end_facet() ;
            } else if(keyword == "f") {
                std::vector<int> cur_facet ;
                while(!in.eol()) {
                    std::string s ;
                    in >> s ;
                    if(s.length() > 0) {
                        std::istringstream v_input(s) ;
                        unsigned int index ;
                        v_input >> index ;
                        if(index < 1 || index > vertex.size()) {
                            std::cerr << "Out of bounds vertex index" 
                                      << std::endl ;
                        } else {
                            cur_facet.push_back(index - 1) ;
                        }
                        char c ;
                        v_input >> c ;
                        // tex vertex, ignored
                        if(c == '/') {
                            v_input >> index ;
                        }
                    }
                }
                if(cur_facet.size() < 3) {
                    std::cerr << "facet with less than 3 vertices, ignoring" 
                              << std::endl ;
                } else {
                    unsigned int f = nb_facets() ;
                    begin_facet() ;
                    for(unsigned int i=0; i<cur_facet.size(); i++) {
                        unsigned int v = cur_facet[i] ;
                        add_vertex(VertexEdge(vertex[v])) ;
                        top_vertex().set_flag(VertexEdge::ORIGINAL) ;
                        vertex_index_.push_back(v) ;
                        star[v].push_back(f) ;
                    }
                    end_facet() ;
                } 
            }
        }
    
        original_vertices_.resize(vertex.size());
        std::copy(vertex.begin(), vertex.end(), original_vertices_.begin());
        
        // Step 2: compute facet adjacencies
        for(unsigned int f=0; f<nb_facets(); f++) {
            unsigned int facet_base = facet_begin(f) ;
            unsigned int facet_n = facet_size(f) ;

            for(unsigned int i=0; i<facet_n; i++) {
                unsigned int v1 = facet_base + i ;
                unsigned int v2 = facet_base + ((i + 1) % facet_n) ;
                unsigned int gv1 = vertex_index_[v1] ;
                unsigned int gv2 = vertex_index_[v2] ;
                const std::vector<unsigned int>& S = star[gv1] ;
                for(unsigned int k=0; k<S.size(); k++) {
                    unsigned int g = S[k] ;
                    if(
                        g != f && has_edge(
                            vertex_index_, 
                            facet_begin(g), facet_end(g), gv2, gv1
                        )
                    ) {
                        this->vertex(v1).f = g ;
                        break ;
                    }
                }
            }
        }

        // Step 3: assign facet ids
        for(unsigned int f=0; f<nb_facets(); f++) {
            facet_info(f).id = f ;
        }

        // Step 4: initialize symbolic information
        init_symbolic_vertices() ;

        // Just for checking
        unsigned int nb_borders = 0 ;
        for(unsigned int i=0; i<nb_vertices(); i++) {
            if(this->vertex(i).f < 0) {
                nb_borders++ ;
            }
        }
        double vol = signed_volume() ;
        orientation_ = (vol > 0.0) ;
        std::cerr << "Mesh loaded, nb_facets = " << nb_facets() 
                  << " nb_borders = " << nb_borders 
                  << " signed volume = " << vol
                  << std::endl ;
        if(!orientation_ && nb_borders == 0) {
            std::cerr << " WARNING ! orientation is negative"
                      << std::endl ;
        }
        return nb_borders ;
    }
    
    void Mesh::init_symbolic_vertices() {
        for(unsigned int f=0; f<nb_facets(); f++) {
            for(unsigned int i1=facet_begin(f); i1<facet_end(f); i1++) {
                unsigned int i2=i1+1 ;
                if(i2 == facet_end(f)) { i2 = facet_begin(f) ; }
                // Note: Here we compute vertex(i2).sym 
                // (and we do not touch vertex(i1).sym).
                vertex(i2).sym.clear() ;
                vertex(i2).sym.add_boundary_facet(f) ;
                if(vertex(i1).f >= 0) {
                    vertex(i2).sym.add_boundary_facet(vertex(i1).f) ;
                } else {
                    // "Virtual" boundary facet: indicates edge (i1,i2 = i1 \oplus 1)
                    vertex(i2).sym.add_boundary_facet(
                        (i1 + nb_facets()) - facet_begin(f)
                    ) ;
                }
                if(vertex(i2).f >= 0) {
                    vertex(i2).sym.add_boundary_facet(vertex(i2).f) ;
                } else {
                    // "Virtual" boundary facet: indicates edge (i2,i2 \oplus 1)
                    vertex(i2).sym.add_boundary_facet(
                        (i2 + nb_facets()) - facet_begin(f) 
                    ) ;
                }
            }
        }        
    }

    void Mesh::save(const std::string& filename) {
        std::ofstream out(filename.c_str()) ;
        for(unsigned int i=0; i<nb_vertices(); i++) {
            out << "v " << vec3(vertex(i)) << std::endl ;
        }
        for(unsigned int f=0; f<nb_facets(); f++) {
            out << "f " ;
            for(unsigned int i=facet_begin(f); i<facet_end(f); i++) {
                out << (i+1) << " " ;
            }
            out << std::endl ;
        }
    }

    class ComputeVolume {
    public:
        ComputeVolume(double& result) : result_(result) {
            result = 0.0 ;
        }
        void operator()(const vec3& p1, const vec3& p2, const vec3& p3) const {
            result_ += mixed_product(p1, p2, p3) / 6.0 ;
        }
    private:
        double& result_ ;
    } ;

    double Mesh::signed_volume() const {
        double result = 0.0 ;
        for_each_triangle(ComputeVolume(result)) ;
        return result ;
    }

    class ComputeArea {
    public:
        ComputeArea(double& result) : result_(result) {
            result = 0.0 ;
        }
        void operator()(const vec3& p1, const vec3& p2, const vec3& p3) const {
            result_ += length(cross(p2-p1,p3-p1)) / 2.0 ;
        }
    private:
        double& result_ ;
    } ;

    double Mesh::area() const {
        double result = 0.0 ;
        for_each_triangle(ComputeArea(result)) ;
        return result ;
    }

//=================================================================================

#ifdef WITH_METIS
    void Mesh::partition(
        unsigned int nb_parts, std::vector<Mesh>& parts
    ) const {
        // Step 0: allocate arrays for METIS
        //   (note: if the mesh has borders, this allocates
        //     slightly more memory than necessary, 
        //     but this is not really a problem since
        //     this is deallocated right after...)
        int n = nb_facets() ;
        int nnz = nb_vertices() ;
        int* xadj = new int[n+1] ;
        int* adjncy = new int[nnz] ;
        int* part = new int[n] ;
        unsigned int* old2new = new unsigned int[n] ;

        // Step 1: copy facet graph and convert into METIS format
        //  (and skip borders !)
        int count = 0 ;
        for(int f=0; f<n; f++) {
            xadj[f] = count ;
            for(unsigned int i=facet_begin(f) ; i<facet_end(f); i++) {
                if(vertex(i).f >= 0) {
                    adjncy[count] = vertex(i).f ;
                    count++ ;
                }
            }
        }
        xadj[n] = count ;

        // Step 2: call METIS
        int options[5] ;
        options[0] = 0 ; // use default values
        int zero = 0;
        int edgecut = 0;
        int nb_parts = nb_parts_in ;

        METIS_PartGraphRecursive(
            &n, xadj, adjncy,  // The matrix
            nil, nil,          // No vertex weight, no edge weight
            &zero,             // No weights
            &zero,             // C-style indexing
            &nb_parts,
            options,
            &edgecut,
            part
        ) ;

        // Step 3: create parts
        parts.resize(nb_parts) ;
        for(unsigned int f=0; f<nb_facets(); f++) {
            int p = part[f] ;
            old2new[f] = parts[p].nb_facets() ;
            parts[p].begin_facet() ;
            for(unsigned int i=facet_begin(f); i<facet_end(f); i++) {
                parts[p].add_vertex(vertex(i)) ;
            }
            parts[p].end_facet() ;
        }

        // Step 4: translate adjacencies
        for(int p=0; p<nb_parts; p++) {
            for(unsigned int v=0; v<parts[p].nb_vertices(); v++) {
                VertexEdge& ve = parts[p].vertex(v) ;
                if(ve.f >= 0) {
                    if(part[ve.f] == p) {
                        ve.f = old2new[ve.f] ;
                    } else {
                        ve.f = -1 ;
                    }
                }
            }
            // Assign facet ids
            for(unsigned int f=0; f<parts[p].nb_facets(); f++) {
                parts[p].facet_info(f).id = f ;
            }
            //   We do not call init_symbolic_vertices() on part[p],
            // since the original symbolic vertices are copied.
            // This means the symbolic vertices have the global 
            // facet IDs (this is what we need).
            //   In contrast, facet_info(f).id is a local index
            // (relative to this part), this is what the RVD algorithm
            // needs.
        }


        // Cleanup
        delete[] xadj ;
        delete[] adjncy ;
        delete[] part ;
        delete[] old2new ;
    }
#else
    void Mesh::partition(
        unsigned int nb_parts, std::vector<Mesh>& parts
    ) const {
        std::cerr << "Cannot partition, need to recompile with METIS support"
                  << std::endl ;
    }

#endif

}
