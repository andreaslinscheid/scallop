/*	This file SelfEnergy.h is part of scallop.
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
 *  Created on: Nov 24, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_SELFENERGY_H_
#define SCALLOP_GW_FLEX_SELFENERGY_H_

#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/gw_flex/GeneralizedSusceptibility.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class SelfEnergy : 	public GreensFunctionOrbital<T>
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;

	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector VbT;

	void set_to_zero();

	void add_electronic_selfenergy(
			GreensFunctionOrbital<T> const& gf,
			GeneralizedSusceptibility<T> const& sf);

	void linear_interpolate_time( VbT const & previousGrid, VbT const & presentGrid );

private:

	void v_matrix(size_t j, size_t a1, size_t s1, size_t a2, size_t s2, T & prefactor) const;

	void transform_2pto_v(size_t j, size_t jp, size_t nOrb,
			MemoryLayout const& twoPtLayout,
			T const * old2PtObj,
			T * transformed2PtObj) const;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/SelfEnergy.hpp"
#endif /* SCALLOP_GW_FLEX_SELFENERGY_H_ */
