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


#ifndef __GEEX_BASICS_LINE_STREAM__
#define __GEEX_BASICS_LINE_STREAM__

#include <iostream>
#include <sstream>

namespace Geex {

	/**
	 * LineInputStream facilitates parsing
	 * ASCII files line by line.
	 */
	class LineInputStream {
	public:
		LineInputStream(std::istream& in) : in_(in), line_in_(0) {   }
		~LineInputStream() { delete line_in_; line_in_ = 0; }
		bool eof() const { return in_.eof(); }
		bool eol() const { return line_in_ == 0 || line_in_->eof(); }
		bool ok() const {
			//return in_ != 0;
			return in_.bad();
		}
		void get_line() {
			in_.getline(buffer_, 65536); delete line_in_; line_in_ = new std::istringstream(buffer_);
		}
		std::istream& line() { return *line_in_; }
		const char* current_line() const { return buffer_; }
		template <class T> LineInputStream& operator>>(T& param) { *line_in_ >> param; return *this; }
	private:
		std::istream& in_;
		std::istringstream* line_in_;
		char buffer_[65536];
	};

}

#endif
