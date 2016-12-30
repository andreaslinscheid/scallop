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
#include "scallop/gw_flex/MemoryLayout.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class GreensFunctionOrbital : public MatsubaraImagTimeFourierTransform<T>,
							  private MemoryLayout
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;

	using MemoryLayout::get_nOrb;

	GreensFunctionOrbital();

	void initialize(
			size_t dimImTime,
			std::vector<size_t> gridDims,
			size_t orbitalDim,
			bool initialInTimeDomain,
			bool initialInReciprocalDomain,
			typename auxillary::TemplateTypedefs<T>::scallop_vector const& data,
			bT chemPot = bT(0));

	void alter_chem_pot_no_shift( bT newChemicalPotential);

	void set_chem_pot( bT newChemicalPotential);

	bT get_chem_pot() const;

	T operator() (
			size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2) const;

	T & operator() (
			size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2);

	T operator() (
			size_t ig, size_t it, size_t m1,  size_t m2) const;

	T & operator() (
			size_t ig, size_t it, size_t m1,  size_t m2);

	void chem_pot_adj_local_part(  bT diffMu, V & localPart) const;

	using MatsubaraImagTimeFourierTransform<T>::transform_itime_Mfreq;

	void transform_itime_Mfreq_subtract(
			bT invT,
			GreensFunctionOrbital<T> const& ksGFTime,
			GreensFunctionOrbital<T> const& ksGFFreq );

private:

	bT chemPot_ = 0;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/GreensFunctionOrbital.hpp"
#endif /* SCALLOP_GW_FLEX_GREENSFUNCTIONORBITAL_H_ */
