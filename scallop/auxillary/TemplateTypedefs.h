/*	This file TemplateTypedefs.h is part of scallop.
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
 *  Created on: Nov 8, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_AUXILLARY_TEMPLATETYPEDEFS_H_
#define SCALLOP_AUXILLARY_TEMPLATETYPEDEFS_H_

#include "scallop/auxillary/AlignmentAllocator.h"
#include <vector>

namespace scallop
{
namespace auxillary
{

template<typename T>
class TemplateTypedefs
{
public:
	typedef std::vector<T, auxillary::AlignmentAllocator<T,32> > scallop_vector;
//	typedef std::vector<T> scallop_vector;

	typedef T (*NambuSpinPauliMatrix)(size_t,size_t,size_t,size_t);
};

} /* namespace auxillary */
} /* namespace scallop */

#endif /* SCALLOP_AUXILLARY_TEMPLATETYPEDEFS_H_ */
