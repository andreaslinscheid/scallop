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

#include "scallop/gw_flex/MemoryLayout.h"
#include "scallop/parallel/GridDistribution.h"
#include "scallop/auxillary/TemplateTypedefs.h"
#include <vector>

namespace scallop
{
namespace gw_flex
{

template<typename T>
class UnitaryWannierKSBands : private MemoryLayout
{
public:
	using MemoryLayout::get_nOrb;

	void initialize_identity(
			std::vector<size_t> spaceGrid,
			size_t numOrbitals);

	void initialize(
			std::vector<size_t> spaceGrid,
			size_t numOrbitals,
			typename auxillary::TemplateTypedefs<T>::scallop_vector data);

	T operator() (size_t ik, size_t m1, size_t m2) const;

	T & operator() (size_t ik, size_t m1, size_t m2);

	T operator() (size_t ik, size_t l, size_t a1, size_t s1, size_t i, size_t a2, size_t s2) const;

	T & operator() (size_t ik, size_t l, size_t a1, size_t s1, size_t i, size_t a2, size_t s2);

	T const * read_phs_grid_ptr(size_t ik, size_t m1, size_t m2 ) const;

	T * write_phs_grid_ptr(size_t ik, size_t m1, size_t m2 );

	T const * read_phs_grid_ptr_block(size_t ik) const;

	T * write_phs_grid_ptr_block(size_t ik);

	parallel::GridDistribution<T> const& get_spaceGrid_proc() const;
private:

	parallel::GridDistribution<T> gdistr_;

	typename auxillary::TemplateTypedefs<T>::scallop_vector data_;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/UnitaryWannierKSBands.hpp"
#endif /* SCALLOP_GW_FLEX_UNITARYWANNIERKSBANDS_H_ */
