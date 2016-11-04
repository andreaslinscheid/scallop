/*	This file UnitaryWannierKSBands.hpp is part of scallop.
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
 *  Created on: Nov 3, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/UnitaryWannierKSBands.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void UnitaryWannierKSBands<T>::initialize_identity(
		std::vector<size_t> spaceGrid,
		size_t numOrbitals)
{
	spaceGrid_ = std::move(spaceGrid);
	numOrb_ = numOrbitals;
	numGridPts_ = 1;
	for ( auto p : spaceGrid_ )
		numGridPts_ *= p;

	data_ = std::vector<T>(16*numGridPts_*numOrb_*numOrb_, T(0) );
	for ( size_t ik = 0 ; ik < numGridPts_ ; ++ik)
		for ( size_t m1 = 0 ; m1 < numOrb_*4 ; ++m1)
		{
			(*this)(ik,m1,m1) = 1.0;
		}
}

template<typename T>
T UnitaryWannierKSBands<T>::operator() (size_t ik, size_t m1, size_t m2) const
{
	return *(data_.begin()+this->memory_layout_combined_notation(ik,m1,m2));
}

template<typename T>
T & UnitaryWannierKSBands<T>::operator() (size_t ik, size_t m1, size_t m2)
{
	return *(data_.begin()+this->memory_layout_combined_notation(ik,m1,m2));
}

template<typename T>
typename std::vector<T>::const_iterator
UnitaryWannierKSBands<T>::get_iterator_at(size_t ik, size_t m1, size_t m2) const
{
	return (data_.begin()+this->memory_layout_combined_notation(ik,m1,m2));
}

template<typename T>
size_t UnitaryWannierKSBands<T>::memory_layout(size_t ik, size_t l, size_t a1, size_t s1, size_t i, size_t a2, size_t s2) const
{
	size_t m1 = (l*2+a1)*2+s1;
	size_t m2 = (i*2+a2)*2+s2;
	return this->memory_layout_combined_notation(ik,m1,m2);
}

template<typename T>
size_t UnitaryWannierKSBands<T>::memory_layout_combined_notation(size_t ik, size_t m1, size_t m2) const
{
	auto nC = numOrb_*4;
	size_t ptr_offset = (ik*nC+m1)*nC+m2;

#ifdef DEBUG_BUILD
//	if ( ptr_offset > ( this->end()-this->begin() ))
//		scallop::error_handling::Error("GreensFunctionOrbital: access out of bounds!");
#endif
	return ptr_offset;
}

template<typename T>
size_t UnitaryWannierKSBands<T>::get_num_orbitals() const
{
	return numOrb_;
}

template<typename T>
size_t UnitaryWannierKSBands<T>::get_num_kpts() const
{
	return numGridPts_;
}

template<typename T>
std::vector<size_t> UnitaryWannierKSBands<T>::get_k_grid() const
{
	return spaceGrid_;
}

} /* namespace gw_flex */
} /* namespace scallop */
