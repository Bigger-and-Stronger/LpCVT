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

// This code has been modified by Canjia Huang <canjia7@gmail.com> on 25-3-1.

#include "LpCVT/others/macro.h"

#include <LpCVT/combinatorics/delaunay_CGAL.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <algorithm>

namespace Geex {

	template <class T> inline vec3g<T>
	triangle_normal_with_area(const vec3g<T>& v1, const vec3g<T>& v2, const vec3g<T>& v3) {
		return cross(v2 - v1, v3 - v1);
	}

	template <class POINT> inline plane3 gx_bisector(
		const POINT& p1, const POINT& p2
	) {
		vec3 gx_p1(p1.cartesian(0), p1.cartesian(1), p1.cartesian(2));
		vec3 gx_p2(p2.cartesian(0), p2.cartesian(1), p2.cartesian(2));
		vec3 N = gx_p1 - gx_p2;
		return Geex::plane3(0.5 * (gx_p1 + gx_p2), N);
	}

	namespace MyCGALStuff {

		// ------------- My stuff to extend CGAL --------------------------

		// Used to mark each individual halfedge of a tetrahedron.
		// Each halfedge is identified by its two vertices's local
		// indices (in the range 0..3). Halfedge flags are packed
		// in a 16 bits integer (there are 12 halfedges in a tet).
		class HalfedgeFlags {
		public:
			HalfedgeFlags() : flags_(0) { }
			void mark(unsigned int i, unsigned int j) {
				flags_ |= mask_[i][j];
			}
			bool is_marked(unsigned int i, unsigned int j) const {
				return ((flags_ & mask_[i][j]) != 0);
			}
			void clear() { flags_ = 0; }
		private:
			unsigned short flags_; // 16 bits
			static unsigned short mask_[4][4];
		};

		// Pre-computed masks for marking and testing 
		// an individual halfedge flag.
		unsigned short HalfedgeFlags::mask_[4][4] = {
			{0,      1,    2,   4},
			{8,      0,   16,  32},
			{64,   128,    0, 256},
			{512, 1024, 2048,   0}
		};

		// ------------- CGAL stuff ---------------------------------        

		// ----------------------- A CGAL::Vertex with decoration ------------------
		template < class Gt, class Vb = CGAL::Triangulation_vertex_base_3<Gt> >
		class Vertex : public  Vb {
			typedef Vb superclass;
		public:
			typedef typename Vb::Vertex_handle      Vertex_handle;
			typedef typename Vb::Cell_handle        Cell_handle;
			typedef typename Vb::Point              Point;

			template < typename TDS2 >
			struct Rebind_TDS {
				typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
				typedef Vertex<Gt, Vb2> Other;
			};

		public:
			Vertex() : superclass(), id(-1) {}
			Vertex(const Point& p) : superclass(p), id(-1) {}
			Vertex(const Point& p, Cell_handle f) : superclass(f, p), id(-1) {}
			Vertex(Cell_handle f) : superclass(f), id(-1) {}

			int id;
		};


		// ----------------------- A CGAL::Cell with decoration ------------------

		template < class Gt, class Cb = CGAL::Delaunay_triangulation_cell_base_3<Gt> >
		class Cell : public Cb {
			typedef Cb superclass;
		public:
			typedef typename Cb::Vertex_handle      Vertex_handle;
			typedef typename Cb::Cell_handle        Cell_handle;
			template < typename TDS2 >
			struct Rebind_TDS {
				typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
				typedef Cell<Gt, Cb2> Other;
			};


			Cell() : superclass() { clear_cicl(); }
			Cell(
				Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
			) : superclass(v0, v1, v2, v3) {
				clear_cicl();
			}

			Cell(
				Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
				Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3
			) : superclass(v0, v1, v2, v3, n0, n1, n2, n3) {
				clear_cicl();
			}

			// Clears Circular Incident Cell List
			void clear_cicl() {
				next_around_vertex_[0] = 0;
				next_around_vertex_[1] = 0;
				next_around_vertex_[2] = 0;
				next_around_vertex_[3] = 0;
			}

			// next_around_vertex stores a circular list
			// of cells incident to a given vertex.
			// It can be traversed as follows:
			// Cell_handle c = v->cell() ;
			// do {
			//   do something with c
			//   c = c->next_around_vertex(c->index(v)) ;
			// } while(c != v->cell()) ;

			Cell_handle next_around_vertex(unsigned int i) const {
				return next_around_vertex_[i];
			}

			void set_next_around_vertex(unsigned int i, Cell_handle c) {
				next_around_vertex_[i] = c;
			}

			HalfedgeFlags& halfedge_flags() { return halfedge_flags_; }


		private:
			Cell_handle next_around_vertex_[4];
			HalfedgeFlags halfedge_flags_;
		public:
			vec3 circumcenter_;
		};

		typedef double Coord_type;
		typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
		typedef Vertex<K> Vb;
		typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb> Vbh;
		//typedef Vb Vbh;
		typedef Cell<K>   Cb;
		typedef CGAL::Triangulation_data_structure_3<Vbh, Cb> TDS;
		typedef CGAL::Delaunay_triangulation_3<K, TDS> TRI;

		typedef K::Point_3    Point;
		typedef K::Vector_3   CGAL_Vector;


		inline vec3 to_geex(const Point& p) {
			return vec3(
				p.cartesian(0), p.cartesian(1), p.cartesian(2)
			);
		}

		class DelaunayBase : public CGAL::Triangulation_hierarchy_3<TRI> {
			// class DelaunayBase : public TRI {
		public:
			std::vector<Vertex_handle> vertices_;

			// Initialize circular incident cell lists.
			void init_cicl() {
				for (Vertex_iterator it = vertices_begin(); it != vertices_end(); it++) {
					Cell_handle c = it->cell();
					c->set_next_around_vertex(c->index(it), c);
				}
				for (Cell_iterator it = cells_begin(); it != cells_end(); it++) {
					for (unsigned int i = 0; i < 4; i++) {
						Vertex_handle v = it->vertex(i);
						if (v->cell() != it) {
							Cell_handle c1 = v->cell();
							unsigned int iv1 = c1->index(v);
							Cell_handle c2 = c1->next_around_vertex(iv1);
							c1->set_next_around_vertex(iv1, it);
							it->set_next_around_vertex(i, c2);
						}
					}
				}
			}

			bool edge_is_marked(const Edge& e) const {
				return e.first->halfedge_flags().is_marked(e.second, e.third);
			}

			// Marks all the halfedges in the tets incident to edge e.
			//   Note1: orientation is taken into account, i.e.
			// (v1,v2) is different from (v2,v1) and can be marked independently.
			//   Note2: we use our next_around_halfedge static table rather than
			// CGAL's circulator since doing so enables us to keep track of local
			// indices and saves some computations (see timings below).
			void mark_edge(const Edge& e) {
				Cell_handle c = e.first;
				unsigned int iv1 = e.second;
				unsigned int iv2 = e.third;
				Vertex_handle v1 = c->vertex(iv1);
				Vertex_handle v2 = c->vertex(iv2);
				Cell_handle cir = c;
				do {
					iv1 = cir->index(v1);
					iv2 = cir->index(v2);
					cir->halfedge_flags().mark(iv1, iv2);
					cir = cir->neighbor(next_around_halfedge[iv1][iv2]);
				} while (cir != c);
			}

		public:
			static unsigned int next_around_halfedge[4][4];
		};

		// Pre-computed table for turning around halfedges in
		// a tetrahedron. Given a cell c and an halfedge (v1, v2), 
		// c->neighbor(next_around_halfedge[v1][v2]) is the cell
		// adjacent to c on the left of (v1,v2).
		// Diagonal entries are not supposed to be accessed.
		unsigned int DelaunayBase::next_around_halfedge[4][4] = {
			{5, 3, 1, 2},
			{2, 5, 3, 0},
			{3, 0, 5, 1},
			{1, 2, 0, 5}
		};


	} // end namespace CGAL 


	class Delaunay_CGAL::Implementation : public ::Geex::MyCGALStuff::DelaunayBase {
	public:
		typedef ::Geex::MyCGALStuff::DelaunayBase baseclass;
	};

	Delaunay_CGAL::Delaunay_CGAL() {
		impl_ = new Implementation;
	}

	Delaunay_CGAL::~Delaunay_CGAL() {
		delete impl_;
	}

	void Delaunay_CGAL::set_vertices(const unsigned int n, const double* xyz) {
		vertices_ = xyz;
		nb_vertices_ = n;
		impl_->clear();

		unsigned int dbl_count = 0;


		if (static_cast<int>(impl_->vertices_.size()) != n)
			impl_->vertices_.resize(n);


		Implementation::Vertex_handle v;
		for (unsigned int i = 0; i < n; i++) {
			v = impl_->insert(Implementation::Point(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]));
			if (v->id != -1) {
				dbl_count++;
				v = nullptr;
			}
			else {
				v->id = i;
			}
			impl_->vertices_[i] = v;
		}

		if (dbl_count > 0) {
			WARNING("Delaunay Warning: encountered " << dbl_count << " duplicated vertices");
		}

		skeleton_dirty_ = true;
	}

	unsigned int Delaunay_CGAL::nearest_vertex_id(double x, double y, double z) const {
		return impl_->nearest_vertex(Implementation::Point(x, y, z))->id;
	}

	// Skeleton timings: 10 Lloyd iterations on pegaso (one vertex per facet)
	// Version 1: 9.12 (uses stack + set to traverse tets, and checks vrtx presence in skel)
	// Version 2: 5.75 (uses incident tet lists and checks vrtx presence in skel)
	// Version 3: 5.00 (uses incident tet lists stored in cells and checks vrtx presence in skel)
	// Version 4: 5.48 (uses incident tet lists stored in cells and marks halfedges)
	// Version 5: 4.82 (same as v.4 but uses next_around_halfedge() instead of CGAL's circulator)

	void Delaunay_CGAL::update_skeleton() const {
		skeleton_dirty_ = false;
		// This is Version 5 (2x faster than Version 1 :-)

		// Compute Circular Incident Cell Lists
		impl_->init_cicl();

		skeleton_.clear();
		for (unsigned int ivg = 0; ivg < impl_->vertices_.size(); ivg++) {
			Implementation::Vertex_handle it = impl_->vertices_[ivg];


			//   This happens when the same vertex was inserted twice
			// in Delaunay (due to two vertices at the same geometric location).
			if (it == 0) { continue; }

			//   Skip nil vertices by inserting empty stars in the
			// Delaunay skeleton (to take into account skipped duplicated
			// vertices, see above).
			while ((int)skeleton_.nb_vertices() < it->id) {
				skeleton_.begin_star();
				// No need to bark, we already said something when we inserted the vertex.
				// std::cerr << "Delaunay Skel Warning: Inserted empty star" << std::endl ;
				skeleton_.end_star();
			}

			skeleton_.begin_star();
			// Traverse incident cell list
			Implementation::Cell_handle c = it->cell();
			do {
				unsigned int ivit = c->index(it);

				// In the current cell, test all edges incident to current vertex 'it'
				for (unsigned int iv = 0; iv < 4; iv++) {
					Implementation::Vertex_handle neigh = c->vertex(iv);

					// Skip (it,it) edge and infinite edges
					if (ivit != iv && !impl_->is_infinite(neigh)) {
						Implementation::Edge e(c, ivit, iv);

						// Skip edges that are already marked
						if (!impl_->edge_is_marked(e)) {

							// Mark current edge
							impl_->mark_edge(e);

							// Add current neighbor and bisector to current star.
							skeleton_.add_to_star(
								neigh->id,
								gx_bisector(it->point(), neigh->point())
							);
						}
					}
				}
				c = c->next_around_vertex(c->index(it));
			} while (c != it->cell());
			skeleton_.end_star();
		}

	}

	void Delaunay_CGAL::get_tetras(std::vector<int>& tetras, bool finite_only) {
		tetras.clear();

		if (finite_only) {
			tetras.reserve(4 * impl_->number_of_finite_cells());
			for (
				Implementation::Finite_cells_iterator it = impl_->finite_cells_begin();
				it != impl_->finite_cells_end(); it++
				) {
				for (unsigned int i = 0; i < 4; i++) {
					tetras.push_back(it->vertex(i)->id);
				}
			}
		}
		else {
			tetras.reserve(4 * impl_->number_of_cells());
			for (
				Implementation::Cell_iterator it = impl_->cells_begin();
				it != impl_->cells_end(); it++
				) {
				for (unsigned int i = 0; i < 4; i++) {
					tetras.push_back(
						impl_->is_infinite(it->vertex(i)) ? -1 : it->vertex(i)->id
					);
				}
			}
		}
	}

	void Delaunay_CGAL::get_facets(std::vector<int>& facets, bool finite_only)
	{
		static int idlist[4][3] = { {1, 2, 3}, {2, 0, 3}, {3, 0, 1}, {0, 2, 1} };


		facets.clear();
		if (finite_only) {
			facets.reserve(3 * impl_->number_of_finite_facets());
			for (
				Implementation::Finite_facets_iterator it = impl_->finite_facets_begin();
				it != impl_->finite_facets_end(); it++
				) {
				Implementation::Cell_handle ch = it->first;
				for (unsigned int i = 0; i < 3; i++) {
					facets.push_back(ch->vertex(idlist[it->second][i])->id);
				}
			}
		}
		else {
			facets.reserve(3 * impl_->number_of_facets());
			for (
				Implementation::Facet_iterator it = impl_->facets_begin();
				it != impl_->facets_end(); it++
				) {
				Implementation::Cell_handle ch = it->first;
				for (unsigned int i = 0; i < 3; i++) {
					Implementation::Vertex_handle vh = ch->vertex(idlist[it->second][i]);
					facets.push_back(
						impl_->is_infinite(vh) ? -1 : vh->id
					);
				}
			}
		}
	}


	void Delaunay_CGAL::get_voronoi_cell(unsigned int iv, VoronoiCell& result, bool geometry) {
		result.clear();

		// update_skeleton creates the circular incident cell lists that we use.
		if (skeleton_dirty_) { update_skeleton(); }

		Implementation::Vertex_handle v = impl_->vertices_[iv];

		assert(v != 0); // Did we have two vertices at same location ?
		// Rem: in this case, we are not likely to ask for the Voronoi cell of the 
		// second one, since it does not generate any intersection.

		// Step 1: clear edge flags 
		Implementation::Cell_handle c = v->cell();
		do {
			c->halfedge_flags().clear();
			c = c->next_around_vertex(c->index(v));
		} while (c != v->cell());

		// Step 2: traverse Delaunay edges
		c = v->cell();
		do {
			unsigned int ivit = c->index(v);
			// In the current cell, test all edges incident to current vertex 'it'
			for (unsigned int iv = 0; iv < 4; iv++) {
				Implementation::Vertex_handle neigh = c->vertex(iv);

				// Skip (it,it) edge and dual facet at infinity
				if (ivit != iv && !impl_->is_infinite(neigh)) {
					Implementation::Edge e(c, ivit, iv);

					// Skip edges that are already marked
					if (!impl_->edge_is_marked(e)) {

						result.begin_facet(neigh->id);
						{
							Implementation::Cell_handle c = e.first;
							unsigned int iv1 = e.second;
							unsigned int iv2 = e.third;
							Implementation::Vertex_handle v1 = c->vertex(iv1);
							Implementation::Vertex_handle v2 = c->vertex(iv2);
							Implementation::Cell_handle cir = c;
							do {
								iv1 = cir->index(v1);
								iv2 = cir->index(v2);
								cir->halfedge_flags().mark(iv1, iv2);
								unsigned int f = Implementation::next_around_halfedge[iv1][iv2];
								unsigned int iv3;
								for (iv3 = 0; iv3 < 4; iv3++) {
									if (iv3 != iv1 && iv3 != iv2 && iv3 != f) {
										break;
									}
								}
								assert(iv3 != 4);
								Implementation::Vertex_handle e_v = cir->vertex(iv3);
								int e_bisector = impl_->is_infinite(e_v) ? -1 : e_v->id;

								bool c_is_infinite = impl_->is_infinite(cir);
								if (geometry) {
									vec3 ccenter(0, 0, 0);
									if (c_is_infinite) {
										for (unsigned int iv = 0; iv < 4; iv++) {
											if (impl_->is_infinite(cir->vertex(iv))) {
												ccenter = 10.0 * triangle_normal_with_area(
													MyCGALStuff::to_geex(cir->vertex((iv + 1) % 4)->point()),
													MyCGALStuff::to_geex(cir->vertex((iv + 2) % 4)->point()),
													MyCGALStuff::to_geex(cir->vertex((iv + 3) % 4)->point())
												);
												if (iv % 2 == 0) {
													ccenter = -ccenter;
												}
												break;
											}
										}
									}
									else {
										ccenter = MyCGALStuff::to_geex(impl_->dual(cir));
									}
									result.add_to_facet(
										e_bisector, ccenter, c_is_infinite
									);
								}
								else {
									result.add_to_facet(
										e_bisector, c_is_infinite
									);
								}
								cir = cir->neighbor(f);
							} while (cir != c);
						}
						result.end_facet();
					}
				}
			}
			c = c->next_around_vertex(c->index(v));
		} while (c != v->cell());
	}


}

CGAL::Delaunay_triangulation_3<Geex::MyCGALStuff::K, Geex::MyCGALStuff::TDS, CGAL::Default, CGAL::Default> x;

CGAL::internal::CC_iterator<CGAL::Compact_container<Geex::MyCGALStuff::Cell<Geex::MyCGALStuff::K, CGAL::Triangulation_cell_base_3<Geex::MyCGALStuff::K, CGAL::Triangulation_ds_cell_base_3<CGAL::Triangulation_data_structure_3<Geex::MyCGALStuff::Vbh, Geex::MyCGALStuff::Cb, CGAL::Sequential_tag>>>>, CGAL::Default, CGAL::Default, CGAL::Default>, false>* xx;
