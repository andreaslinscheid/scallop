/*	This file GreensFunctionOrbital.hpp is part of scallop.
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
 *  Created on: Nov 1, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/error_handling/error_handling.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
GreensFunctionOrbital<T>::GreensFunctionOrbital()
	: MatsubaraImagTimeFourierTransform<T>(true)
{

}

template<typename T>
void GreensFunctionOrbital<T>::initialize(
		size_t dimImTime,
		std::vector<size_t> gridDims,
		size_t orbitalDim,
		bool initialInFreqDomain,
		bool initialInReciprocalDomain,
		std::vector<T> data)
{
	orbitalDim_ = orbitalDim;

	MatsubaraImagTimeFourierTransform<T>::initialize(
			dimImTime,
			gridDims,
			orbitalDim*orbitalDim*16,
			initialInFreqDomain,
			initialInReciprocalDomain,
			std::move(data));
}

template<typename T>
T GreensFunctionOrbital<T>::operator() (
		size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2) const
{
	return *(this->data_begin_modify() + this->memory_layout(ig,it,l1,a1,s1,l2,a2,s2) );
}

template<typename T>
T & GreensFunctionOrbital<T>::operator() (
		size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2)
{
	return *(this->data_begin_modify() + this->memory_layout(ig,it,l1,a1,s1,l2,a2,s2) );
}

template<typename T>
size_t GreensFunctionOrbital<T>::memory_layout(
		size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2) const
{
	auto nM = this->get_num_time();
	auto nO = orbitalDim_;
	size_t ptr_offset = ((((((ig*nM+it)*nO+l1)*2+a1)*2+s1)*nO+l2)*2+a2)*2+s2;

#ifdef DEBUG_BUILD
	if ( ptr_offset > ( this->end()-this->begin() ))
		scallop::error_handling::Error("GreensFunctionOrbital: access out of bounds!");
#endif
	return ptr_offset;
}

} /* namespace gw_flex */
} /* namespace scallop */
