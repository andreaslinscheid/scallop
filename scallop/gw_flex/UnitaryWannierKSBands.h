/*	This file UnitaryWannierKSBands.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_UNITARYWANNIERKSBANDS_H_
#define SCALLOP_GW_FLEX_UNITARYWANNIERKSBANDS_H_

#include <vector>

namespace scallop
{
namespace gw_flex
{

template<typename T>
class UnitaryWannierKSBands
{
public:

	void initialize_identity(
			std::vector<size_t> spaceGrid,
			size_t numOrbitals);

	size_t get_num_orbitals() const;

	size_t get_num_kpts() const;

	std::vector<size_t> get_k_grid() const;

	T operator() (size_t ik, size_t m1, size_t m2) const;

	T & operator() (size_t ik, size_t m1, size_t m2);

	typename std::vector<T>::const_iterator
		get_iterator_at(size_t ik, size_t m1, size_t m2) const;

	T operator() (size_t ik, size_t l, size_t a1, size_t s1, size_t i, size_t a2, size_t s2) const;

	T & operator() (size_t ik, size_t l, size_t a1, size_t s1, size_t i, size_t a2, size_t s2);

	typename std::vector<T>::const_iterator
		get_iterator_at(size_t ik, size_t l, size_t a1, size_t s1, size_t i, size_t a2, size_t s2) const;
private:

	size_t numOrb_ = 0;

	size_t numGridPts_ = 0;

	std::vector<size_t> spaceGrid_;

	std::vector<T> data_;

	size_t memory_layout(size_t ik, size_t l, size_t a1, size_t s1, size_t i, size_t a2, size_t s2) const;

	size_t memory_layout_combined_notation(size_t ik, size_t m, size_t m2) const;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/UnitaryWannierKSBands.hpp"
#endif /* SCALLOP_GW_FLEX_UNITARYWANNIERKSBANDS_H_ */
