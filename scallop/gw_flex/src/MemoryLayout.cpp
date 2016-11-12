/*	This file MemoryLayout.cpp is part of scallop.
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
 *  Created on: Nov 9, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/MemoryLayout.h"

namespace scallop
{
namespace gw_flex
{

void MemoryLayout::initialize_layout_2pt_obj( size_t numOrbitals )
{
	numOrbitals_ = numOrbitals;
}

size_t MemoryLayout::get_nOrb() const
{
	return numOrbitals_;
}

size_t MemoryLayout::memory_layout_2pt_obj(size_t l1, size_t a1, size_t s1, size_t l2, size_t a2, size_t s2) const
{
	size_t m1 = (l1*2+a1)*2+s1;
	size_t m2 = (l2*2+a2)*2+s2;
	return this->memory_layout_combined_notation_2pt_obj(m1,m2);
}

size_t MemoryLayout::memory_layout_combined_notation_2pt_obj(size_t m1, size_t m2) const
{
	auto nC = numOrbitals_*4;
	return m1*nC+m2;
}

} /* namespace gw_flex */
} /* namespace scallop */
