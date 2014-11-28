/*	This file NambuOffDiagonalEnergyIntegralN.h is part of scallop.
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
 *  Created on: Nov 26, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_ELIASHBERG_NAMBUOFFDIAGONALENERGYINTEGRALN_H_
#define SCALLOP_ELIASHBERG_NAMBUOFFDIAGONALENERGYINTEGRALN_H_

#include "scallop/eliashberg/MatzubaraSplittingVector.h"

namespace scallop {
namespace eliashberg {

template<typename T>
class NambuOffDiagonalEnergyIntegralN : public MatzubaraSplittingVector<T> {

	template<class gapSquareRoot, class FrequencyRenorm, class AsymetricEnergyRenorm>
	void compute(
			gapSquareRoot const& gapSqrtUp,
			gapSquareRoot const& gapSqrtDown,
			FrequencyRenorm const& Z,
			AsymetricEnergyRenorm const& A);
};

} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/src/NambuOffDiagonalEnergyIntegralN.hpp"
#endif /* SCALLOP_ELIASHBERG_NAMBUOFFDIAGONALENERGYINTEGRALN_H_ */
