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

#ifndef __GEEX_CVT_CONTAINERS__
#define __GEEX_CVT_CONTAINERS__

#include <set>
#include <vector>
#include <algorithm>
#include <iostream>

namespace Geex {

    /**
     * Similar to std::set, but with fixed maximum size
     * (and no dynamic memory allocation). Used to store
     * vertices equations (represented as plane indices 
     * triplets).
     */
    template <class T, int DIM> class small_set {
    public:
        typedef small_set<T,DIM> thisclass ;
        typedef T* iterator ;
        typedef const T* const_iterator ;
        typedef T& reference ;
        typedef T value_type ;
        
        small_set() : end_(data_) { }

        small_set(const thisclass& rhs) {
            copy(rhs) ;
        }
        
        thisclass& operator=(const thisclass& rhs) {
            copy(rhs) ; return *this ;
        }

        template <int DIM2> small_set(const small_set<T,DIM2>& rhs) {
            copy(rhs) ;
        }
        
        template <int DIM2> thisclass& operator=(const small_set<T,DIM2>& rhs) {
            copy(rhs) ; return *this ;
        }

        unsigned int size() const { return (unsigned int)(end_ - data_) ; }
        unsigned int capacity() const { return (unsigned int)DIM ; }

        iterator begin() { return data_ ; }
        iterator end() { return end_ ; }
        iterator end_of_storage() { return data_ + DIM ; }

        iterator* end_ptr() { return &end_ ; }

        const_iterator begin() const { return data_ ; }
        const_iterator end() const { return end_ ; }
        const_iterator end_of_storage() const { return data_ + DIM ; }

        iterator insert(const T& x) { 
            return insert(x, find_i(x)) ; 
        }

        iterator insert(const T& x, iterator where) {
            if(where == end()) { *where = x ; grow() ; return where ; }
            if(*where == x) { return where ; }
            grow() ;
            if(where == end() - 1) {
                *where = x ; return where ; 
            }
            for(iterator i = end()-1; i != where; i--) {
                *i = *(i-1) ;
            }
            *where = x ;
            return where ;
        }

        void clear() { end_ = data_ ; }

        iterator find(const T& x) {
            iterator result = find_i(x) ;
            if(*result != x) { result = end() ; }
            return result ;
        }

        const_iterator find(const T& x) const {
            const_iterator result = find_i(x) ;
            if(*result != x) { result = end() ; }
            return result ;
        }

        void push_back(const T& x) {
            *end_ = x ;
            grow() ;
        }

        void print(std::ostream& out) const {
            out << "[ " ;
            for(const_iterator it=begin(); it!=end(); it++) {
                out << *it << " " ;
            }
            out << "]" ;
        }

        T& operator[](int i) {
            return begin()[i] ;
        }

        const T& operator[](int i) const {
            return begin()[i] ;
        }

    protected:

        void grow() {
            end_++ ;
        }


        template <int DIM2> void copy(const small_set<T,DIM2>& rhs) {
            end_ = data_ ;
            for(typename small_set<T,DIM2>::const_iterator 
                    it=rhs.begin(); it!=rhs.end(); it++) {
                push_back(*it) ;
            }
        }


        // Note: maybe we should start from end() instead of begin()
        // since negative indices are inserted first.
        iterator find_i(const T& x) {
            iterator result = begin() ;
            while(result != end() && *result < x) {
                result++ ;
            }
            return result ;
        }

        const_iterator find_i(const T& x) const {
            const_iterator result = begin() ;
            while(result != end() && *result < x) {
                result++ ;
            }
            return result ;
        }

    protected:
        T data_[DIM] ;
        iterator end_ ;
    } ;

    template <class T, int DIM> inline std::ostream& operator<<(
        std::ostream& out,
        const small_set<T, DIM>& S
    ) {
        S.print(out) ; return out ;
    }


    template <class T, int DIM1, int DIM2, int DIM3> inline void sets_intersect(
        const small_set<T,DIM1>& S1, const small_set<T,DIM2>& S2, small_set<T,DIM3>& I
    ) {
        I.clear() ;
        typename small_set<T,DIM1>::const_iterator i1 = S1.begin() ;
        typename small_set<T,DIM2>::const_iterator i2 = S2.begin() ; 
        while(i1 < S1.end() && i2 < S2.end()) {
            if(*i1 < *i2) 
                ++i1 ;
            else if (*i2 < *i1) 
                ++i2 ;
            else {
                I.push_back(*i1) ;
                ++i1 ;
                ++i2 ;
            }
        }
    }

    //--------------------------------------------------------------------------------

    template <class T, int DIM> inline std::ostream& operator<<(
        std::ostream& out,
        const std::set<T>& S
    ) {
        out << "[ " ;
        for(typename std::set<T>::const_iterator it = S.begin(); it != S.end(); it++) {
            out << *it << " " ;
        }
        out << "]" ;
        return out ;
    }
    
    template <class T> inline void sets_intersect(
        const std::set<T>& S1, const std::set<T>& S2, std::set<T>& I
    ) {
        I.clear() ;
        std::set_intersection(
            S1.begin(), S1.end(), S2.begin(), S2.end(), std::inserter(I, I.begin())
        ) ;
    }
    

    //----------------------------------------------------------------------------
}


#endif
