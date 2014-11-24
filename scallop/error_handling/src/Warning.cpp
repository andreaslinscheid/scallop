/*	This file Warning.cpp is part of scallop.
 *
 *  scallop is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  scallop is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with scallop.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Nov 24, 2014
 *      Author: Andreas Linscheid
 */

#include "scallop/error_handling/Warning.h"

namespace scallop {
namespace error_handling {

Warning::Warning() : _buffer(), _sstrBuff(_buffer) { };

void Warning::print() const {
	if ( _sstrBuff.str().empty() )
		return;
	std::cerr << "=====================================================\n"
				 "||WARNING : ||\n" <<  _sstrBuff.str() <<
				 "=====================================================\n" <<std::endl;
}

Warning::~Warning() {
	this->print();
}

Warning::Warning(std::string const& message) :
		_buffer(message), _sstrBuff(_buffer)  { };

} /* namespace error_handling */
} /* namespace scallop */
