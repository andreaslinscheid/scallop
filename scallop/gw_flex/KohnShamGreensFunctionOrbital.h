/*	This file KohnShamGreensFunctionOrbital.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_KOHNSHAMGREENSFUNCTIONORBITAL_H_
#define SCALLOP_GW_FLEX_KOHNSHAMGREENSFUNCTIONORBITAL_H_

#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/gw_flex/UnitaryWannierKSBands.h"
#include "scallop/auxillary/TypeMapComplex.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class KohnShamGreensFunctionOrbital : public GreensFunctionOrbital<T>
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	void set_in_time_space(
			UnitaryWannierKSBands<T> const& unitaryWannierBands,
			std::vector<bT> const& KSBands,
			size_t timeDim,
			bT invTemp);

	void set_in_frequency_space(
			UnitaryWannierKSBands<T> const& unitaryWannierBands,
			std::vector<bT> const& KSBands,
			size_t freqDim,
			bT invTemp);
private:

	void set_in_both_spaces(
			UnitaryWannierKSBands<T> const& unitaryWannierBands,
			std::vector<bT> const& KSBands,
			size_t timeDim,
			bT invTemp,
			bool timeSpace);
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/KohnShamGreensFunctionOrbital.hpp"
#endif /* SCALLOP_GW_FLEX_KOHNSHAMGREENSFUNCTIONORBITAL_H_ */
