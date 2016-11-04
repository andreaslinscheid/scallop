/*	This file GreensFunctionOrbital.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_GREENSFUNCTIONORBITAL_H_
#define SCALLOP_GW_FLEX_GREENSFUNCTIONORBITAL_H_

#include "scallop/gw_flex/MatsubaraImagTimeFourierTransform.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class GreensFunctionOrbital : public MatsubaraImagTimeFourierTransform<T>
{
public:

	GreensFunctionOrbital();

	void initialize(
			size_t dimImTime,
			std::vector<size_t> gridDims,
			size_t orbitalDim,
			bool initialInFreqDomain,
			bool initialInReciprocalDomain,
			std::vector<T> data);

	T operator() (
			size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2) const;

	T & operator() (
			size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2);

	T operator() (
			size_t ig, size_t it, size_t m1,  size_t m2) const;

	T & operator() (
			size_t ig, size_t it, size_t m1,  size_t m2);

	typename std::vector<T>::iterator
	get_iterator_at(
			size_t ig, size_t it, size_t m1, size_t m2);

private:

	size_t orbitalDim_ = 0;

	size_t memory_layout(
			size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2) const;

	size_t memory_layout_combined_notation(
			size_t ig, size_t it, size_t m1, size_t m2) const;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/GreensFunctionOrbital.hpp"
#endif /* SCALLOP_GW_FLEX_GREENSFUNCTIONORBITAL_H_ */
