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

template<typename T>
Warning& Warning::operator<< (T const& message) {
	 _sstrBuff << message;
	return *this;
}

} /* namespace error_handling */
} /* namespace scallop */
