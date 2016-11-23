/*	This file MemoryLayout.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_SRC_MEMORYLAYOUT_H_
#define SCALLOP_GW_FLEX_SRC_MEMORYLAYOUT_H_

#include <cstdlib>

namespace scallop
{
namespace gw_flex
{

/**
 * This object defines the local memory layout within a black of data at a point of space and time.
 */
class MemoryLayout
{
public:

	void initialize_layout_2pt_obj( size_t numOrbitals );

	void initialize_layout_4pt_scalar_obj( size_t numOrbitals, size_t spinChargeChannels = 4 );

	void initialize_layout_phonon_prop( size_t numModes );

	size_t memory_layout_2pt_obj(size_t l1, size_t a1, size_t s1, size_t l2, size_t a2, size_t s2) const;

	size_t memory_layout_2pt_obj_nsc(size_t l1, size_t as1, size_t l2, size_t as2) const;

	size_t memory_layout_combined_notation_2pt_obj(size_t m1, size_t m2) const;

	size_t memory_layout_4pt_scalar_obj(size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4) const;

	size_t memory_layout_combined_notation_4pt_scalar_obj(size_t j, size_t jp, size_t m1, size_t m2) const;

	size_t memory_layout_phonon_prop(size_t nu, size_t nup) const;

	size_t get_nOrb() const;

	size_t get_nChnls() const;
private:

	size_t numOrbitals_ = 0;

	size_t spinChargeChannels_ = 0;
};

} /* namespace gw_flex */
} /* namespace scallop */

#endif /* SCALLOP_GW_FLEX_SRC_MEMORYLAYOUT_H_ */
