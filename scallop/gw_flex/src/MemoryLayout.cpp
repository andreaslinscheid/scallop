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
#include <assert.h>

namespace scallop
{
namespace gw_flex
{

void MemoryLayout::initialize_layout_4pt_scalar_obj( size_t numOrbitals, size_t spinChargeChannels )
{
	numOrbitals_ = numOrbitals;
	spinChargeChannels_ = spinChargeChannels;
}

void MemoryLayout::initialize_layout_2pt_obj( size_t numOrbitals )
{
	numOrbitals_ = numOrbitals;
}

void MemoryLayout::initialize_layout_phonon_prop( size_t numModes )
{
	numOrbitals_ = numModes;
}

size_t MemoryLayout::get_nOrb() const
{
	return numOrbitals_;
}

size_t MemoryLayout::get_nChnls() const
{
	return spinChargeChannels_;
}

size_t MemoryLayout::memory_layout_2pt_diagonal(size_t l1, size_t a1, size_t s1) const
{
	size_t as1 = a1*2+s1;
	assert( (a1 < 2) && ( s1 < 2) && (l1 < numOrbitals_) );
	return this->memory_layout_2pt_diagonal_nsc(l1,as1);
}

size_t MemoryLayout::memory_layout_2pt_diagonal_nsc(size_t l1, size_t as1) const
{
	assert( (as1 < 4) && (l1 < numOrbitals_) );
	return as1*numOrbitals_+l1;
}

size_t MemoryLayout::memory_layout_2pt_obj(size_t l1, size_t a1, size_t s1, size_t l2, size_t a2, size_t s2) const
{
	size_t as1 = a1*2+s1;
	size_t as2 = a2*2+s2;
	assert( (a1 < 2) && ( s1 < 2) && (l1 < numOrbitals_) && (a2 < 2) && ( s2 < 2) && (l2 < numOrbitals_));
	return this->memory_layout_2pt_obj_nsc(l1,as1,l2,as2);
}

size_t MemoryLayout::memory_layout_2pt_obj_nsc(size_t l1, size_t as1, size_t l2, size_t as2) const
{
	size_t m1 = as1*numOrbitals_+l1;
	size_t m2 = as2*numOrbitals_+l2;
	assert( (as1 < 4) && (l1 < numOrbitals_) && (as2 < 4) && (l2 < numOrbitals_));
	return this->memory_layout_combined_notation_2pt_obj(m1,m2);
}

size_t MemoryLayout::memory_layout_combined_notation_2pt_obj(size_t m1, size_t m2) const
{
	auto nC = numOrbitals_*4;
	assert( (m1 < nC) && ( m2 < nC) );
	return m1*nC+m2;
}

size_t MemoryLayout::memory_layout_4pt_scalar_obj(size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4) const
{
	size_t m1 = l1*numOrbitals_+l2;
	size_t m2 = l3*numOrbitals_+l4;
	assert( (j < spinChargeChannels_) && ( jp < spinChargeChannels_)
			&& (l1 < numOrbitals_) && (l2 < numOrbitals_) && (l3 < numOrbitals_) && (l4 < numOrbitals_));
	return this->memory_layout_combined_notation_4pt_scalar_obj(j,jp,m1,m2);
}

size_t MemoryLayout::memory_layout_combined_notation_4pt_scalar_obj(size_t j, size_t jp, size_t m1, size_t m2) const
{
	auto nC = numOrbitals_*numOrbitals_;
	assert( (j < spinChargeChannels_) && ( jp < spinChargeChannels_)
			&& (m1 < nC) && (m2 < nC) );
	return ((j*spinChargeChannels_+jp)*nC+m1)*nC+m2;
}

size_t MemoryLayout::memory_layout_phonon_prop(size_t nu, size_t nup) const
{
	assert( (nu < numOrbitals_) && (nup < numOrbitals_) );
	return nu*numOrbitals_ + nup;
}

} /* namespace gw_flex */
} /* namespace scallop */
