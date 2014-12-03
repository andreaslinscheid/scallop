/*	This file GapSquareRoot.h is part of scallop.
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
 *  Created on: Nov 27, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_ELIASHBERG_GAPSQUAREROOT_H_
#define SCALLOP_ELIASHBERG_GAPSQUAREROOT_H_

#include "scallop/eliashberg/MatzubaraSplittingVector.h"
#include "scallop/eliashberg/FrequencyCorrectionZ.h"
#include "scallop/eliashberg/EliashbergGapFunction.h"
#include "scallop/eliashberg/ASymmetricEnergyCorrection.h"

namespace scallop {
namespace eliashberg {

template<typename T, int spin>
class GapSquareRoot : public MatzubaraSplittingVector<T> {

	static_assert( abs(spin) == 1, "Must be specialized with valid spin of +-1" );

public :

	void compute(FrequencyCorrectionZ<T> const& Z,
			EliashbergGapFunction<T> const& Delta,
			ASymmetricEnergyCorrection<T> const& A,
			T inverseTemperature);
};

} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/src/GapSquareRoot.hpp"
#endif /* SCALLOP_ELIASHBERG_GAPSQUAREROOT_H_ */
